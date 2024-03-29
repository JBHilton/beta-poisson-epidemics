'''This file contains all of the datasets we use in our analysis.'''

plague_data = [0]*16+[1]*10+[2]*7+[3]*2+[4]*3+[5]+[6]

mpox_data_g1 = [0]*114+[1]*23+[2]*8+[3]*1+[5]*1 # First generation cases
mpox_data_g2 = [0]*38+[1]*7+[2]*1+[3]*1 # Second generation cases
mpox_data_g3 = [0]*9+[1]+[2] # Third gen cases
mpox_data_g4 = [0]*2+[1] # Fourth gen cases
mpox_data = mpox_data_g1 + \
                 mpox_data_g2 + \
                 mpox_data_g3 + \
                 mpox_data_g4

nigeria_ebola_data = [0]*15+[1]*2+[2]+[3]+[12]

guinea_ebola_data = [1,2,2,5,14,1,4,4,1,3,3,8,2,1,1,4,9,9,1,1,17,2,1,1,1,4,3,3,4,2,
                5,1,2,2,1,9,1,3,1,2,1,1,2]
guinea_ebola_data = guinea_ebola_data + [0] * (152 - len(guinea_ebola_data))

singapore_sars_data = [0]*162+[1]*19+[2]*8+[3]*7+[7]+[12]+[21]+[23]+[40]

sk_mers_data = [38,3,2,1,6,81,2,23,2,1,1,1,5,1,1,1,2,1,1,1]
sk_mers_data = sk_mers_data + [0] * (166 - len(sk_mers_data))

sa_mers_data = [0]*13+[1]*5+[2]*4+[3]+[7]

noro_data = [0]*22+[1]*13+[2]*6+[3]*3+[4]+[5]
