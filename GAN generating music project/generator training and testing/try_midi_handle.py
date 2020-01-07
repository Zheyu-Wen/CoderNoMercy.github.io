
import midi
from matplotlib.pyplot import *

lowerBound = 24
upperBound = 102
span = upperBound - lowerBound


def midi_information(midifile, label, ti, pit, squash=True, span=span):
    pattern = midi.read_midifile(midifile)
    # print(pattern)
    for i in range(len(pattern)):
        one_ti = []
        one_label = []
        one_vel = []
        one_pit = []
        for j in range(len(pattern[i])):
            evt = pattern[i][j]
            if isinstance(evt, midi.NoteEvent):
                if isinstance(evt, midi.NoteOffEvent):
                    one_label.append(0)
                else:
                    one_label.append(1)
                one_ti.append(evt.tick)
                one_pit.append(evt.pitch)
                # one_vel.append(evt.velocity)
        one_label = np.array(one_label)
        #one_vel = np.array(one_vel)
        one_ti = np.array(one_ti)
        one_pit = np.array(one_pit)
        one_label = np.transpose(one_label)
        one_ti = np.transpose(one_ti)
        one_pit = np.transpose(one_pit)
        #one_vel = np.transpose(one_vel)
        label.append(one_label)
        pit.append(one_pit)
        ti.append(one_ti)
        # vel.append(one_vel)

    return label, ti, pit
