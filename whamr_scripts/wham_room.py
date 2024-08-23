import numpy as np
import pyroomacoustics as pra
from pyroomacoustics.parameters import constants
from scipy.signal import resample_poly

class WhamRoom(pra.room.ShoeBox):

    def __init__(self, p, mics, s1, s2, T60, fs=16000,
                 t0=0., sigma2_awgn=None):

        self.T60 = T60
        self.max_rir_len = np.ceil(T60*fs).astype(int)

        volume = p[0]*p[1]*p[2]
        surface_area = 2*(p[0]*p[1] + p[0]*p[2] + p[1]*p[2])
        absorption = 24 * volume * np.log(10.0) / (constants.get('c') * surface_area * T60)

        # minimum max order to guarantee complete filter of length T60
        max_order = np.ceil(T60 * constants.get('c') / min(p)).astype(int)

        super().__init__(p, fs=fs, t0=t0, absorption=absorption,
                         max_order=max_order, sigma2_awgn=sigma2_awgn,
                         sources=None, mics=None)

        self.add_source(s1)
        self.add_source(s2)

        self.add_microphone_array(pra.MicrophoneArray(np.array(mics).T, fs))

    def add_audio(self, s1, s2):
        self.sources[0].add_signal(s1)
        self.sources[1].add_signal(s2)

    def compute_rir(self):

        self.rir = []
        self.rir_early = []
        self.visibility = None

        self.image_source_model()

        for m, mic in enumerate(self.mic_array.R.T):
            h = []
            h_early = []
            for s, source in enumerate(self.sources):
                rir_result = source.get_rir(mic, self.visibility[s][m], self.fs, self.t0)[:self.max_rir_len]
                rir_early = rir_result
                idx = np.argmax(rir_result)
                if self.fs == 16000:
                    rir_early[idx+41:] = 0
                elif self.fs == 8000:
                    rir_early[idx+21:] = 0
                h.append(rir_result)
                h_early.append(rir_early)
            self.rir.append(h)
            self.rir_early.append(h_early)

    def generate_rirs(self):

        original_max_order = self.max_order
        self.max_order = 0

        self.compute_rir()

        self.rir_anechoic = self.rir

        self.max_order = original_max_order

        self.compute_rir()

        self.rir_reverberant = self.rir

    def generate_audio(self, anechoic=False, fs=16000):

        if not self.rir:
            self.generate_rirs()
        if anechoic:
            self.rir = self.rir_anechoic
        else:
            self.rir = self.rir_reverberant
        audio_array = self.simulate(return_premix=True, recompute_rir=False)

        if type(fs) is not list:
            fs_array = [fs]
        else:
            fs_array = fs
        audio_out = []
        for elem in fs_array:
            if type(elem) is str:
                elem = int(elem.replace('k','000'))
            if elem != self.fs:
                assert(self.fs % elem == 0)
                audio_out.append(resample_poly(audio_array, elem, self.fs, axis=2))
            else:
                audio_out.append(audio_array)
        if type(fs) is not list:
            return audio_out[0] # array of shape (n_sources, n_mics, n_samples)
        else:
            return audio_out
        
    def generate_audio_early(self, anechoic=False, fs=16000):

        if not self.rir:
            self.generate_rirs()
        if anechoic:
            self.rir = self.rir_anechoic
        else:
            self.rir = self.rir_early
        audio_array = self.simulate(return_premix=True, recompute_rir=False)

        if type(fs) is not list:
            fs_array = [fs]
        else:
            fs_array = fs
        audio_out = []
        for elem in fs_array:
            if type(elem) is str:
                elem = int(elem.replace('k','000'))
            if elem != self.fs:
                assert(self.fs % elem == 0)
                audio_out.append(resample_poly(audio_array, elem, self.fs, axis=2))
            else:
                audio_out.append(audio_array)
        if type(fs) is not list:
            return audio_out[0] # array of shape (n_sources, n_mics, n_samples)
        else:
            return audio_out