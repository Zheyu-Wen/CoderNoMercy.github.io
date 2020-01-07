# coding = utf-8
# encoding = utf-8
import os
import try_midi_handle as MH
import random
from try_gan_midi import model
import try_note_to_midi as NM
import numpy as np
import midi_state

batch_size = 400
one_out_size = 400
times = 1
lasting = 1
scale = 10

file_name = []
for root, _, _ in os.walk('/home/yht/midi/data'):
    file_name.append(root)

for k in range(len(file_name)):
    for y in range(5):
        label = []
        ti = []
        pit = []
        vel = []
        if k == len(file_name) - 1:
            break
        midifile = file_name[k + 1]
        for _, _, files in os.walk('{}'.format(midifile)):
            for i in range(len(files)):
                [label, ti, pit] = MH.midi_information('{}/{}'.format(midifile, files[i]), label, ti, pit,
                                                       span=MH.span)
            each_all = []
            new_pit = []
            new_label = []
            new_ti = []
            new_vel = []
            new_state = []

            posi = 0
            label = np.array(label)
            ti = np.array(ti)
            pit = np.array(pit)
            length = label.shape[0]
            posi = random.sample([x for x in range(length)], 1)[0]
            print length, posi
            one_git = model(batch_size, one_out_size, times, posi, pit)
            while one_git == 0:
                posi = random.sample([x for x in range(length)], 1)[0]
                one_git = model(batch_size, one_out_size, times, posi, pit)

            length = len(one_git)
            new_pit.append(one_git)
            new_label.append(label[posi][0:length])
            new_state.append(midi_state.count_state(label[posi]))
            new_ti.append(ti[posi][0:length])
            # m,n=new_pit.shape
            # new_new_pit = np.fft.fft(new_pit)
            # plot(np.array(list(range(0,n*m))),new_new_pit)
            # show()
            matrix = [new_pit, new_label, new_ti, new_state]

            NM.note_to_midi(matrix=matrix, name1=k, name2=y)
