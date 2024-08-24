import os

test_mix_scp = 'albert_mix.scp'

test_mix = '/home/nas3/user/albert/graduation/DB/output_v1/mixed/'



tt_mix = open(test_mix_scp,'w')

for i in range(1, 101):
    tt_mix.write(str(i) + ".wav "+ " " +test_mix+str(i)+"/x.wav")
    tt_mix.write('\n')

# for root, dirs, files in os.walk(test_mix):
#     files.sort()
#     for file in files:
#         tt_mix.write(file+" "+root+'/'+file)
#         tt_mix.write('\n')
