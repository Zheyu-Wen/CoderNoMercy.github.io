# coding = utf-8
# encoding = utf-8
from __future__ import division
import midi
from matplotlib.pyplot import *
import random
import os
import pdb
music_combine = {0: [[48, 52, 55], [48, 53, 57]], 1: [[49, 54, 57], [49, 55, 59]], 2: [[50, 54, 57], [50, 55, 59]],
                 3: [[51, 55, 58], [51, 56, 48]], 4: [[52, 57, 49], [52, 57, 49]], 5: [[53, 57, 48], [53, 59, 51]],
                 6: [[54, 58, 94], [54, 59, 51]], 7: [[55, 59, 50], [55, 48, 52]], 8: [[56, 59, 51], [56, 49, 52]],
                 9: [[57, 49, 52], [57, 50, 54]], 10: [[58, 49, 53], [58, 51, 54]], 11: [[59, 51, 54], [59, 52, 56]]}


batch_size = 400
one_out_size = 400
times = 1
lasting_on = 120
lasting_off = 120
lasting = lasting_off + lasting_on


def get_track2(track, pit, music_combine):
    pit = pit % 12
    posi = random.sample([0, 1], 1)[0]
    vel = 60
    print pit, posi
    new_pit = music_combine[pit][posi]
    track.append(midi.NoteOnEvent(tick=240, channel=1, velocity=vel, pitch=new_pit[0]))
    track.append(midi.NoteOnEvent(tick=0, channel=1, velocity=vel, pitch=new_pit[1]))
    track.append(midi.NoteOnEvent(tick=0, channel=1, velocity=vel, pitch=new_pit[2]))
    track.append(midi.NoteOffEvent(tick=240, channel=1, velocity=vel, pitch=new_pit[0]))
    track.append(midi.NoteOffEvent(tick=0, channel=1, velocity=vel, pitch=new_pit[1]))
    track.append(midi.NoteOffEvent(tick=0, channel=1, velocity=vel, pitch=new_pit[2]))
    return track


def find_track(random_bat, track_num):
    training_track = []
    for i in range(4):
        if random_bat[0] < track_num:
            training_track.append(random_bat[0])
            random_bat[0] = random_bat[0] + 1
    return training_track


def note_to_midi(matrix, name1, name2):
    pattern = midi.Pattern()
    handled_pitch = 0
    track_num = len(matrix[0])
    time = 1
    lasttime = 0
    for i in range(track_num):
        if (i == 0):
            vel = 100
        else:
            vel = 40
        count_state = 0
        track = midi.Track()
        track1 = midi.Track()
        pattern.append(track)
        pattern.append(track1)
        pit = matrix[0][i]
        # pdb.set_trace()
        label = matrix[1][i]
        tick = matrix[2][i]

        state = matrix[3][i]
        length = np.shape(pit)[0]
        if len(state) == 0:
            for j in range(length):
                handled_pitch = handled_pitch + 1
                if((handled_pitch % (960 / lasting)) == 0):
                    if handled_pitch < length:
                        track1 = get_track2(track1, pit[handled_pitch], music_combine)
                one_label = label[j]
                one_pit = int(pit[j])
                # one_vel = vel[j]
                one_tick = tick[j]
                if one_label == 0:
                    track.append(midi.NoteOffEvent(tick=(time - lasttime) * lasting_off, channel=i, pitch=one_pit))
                else:
                    track.append(
                        midi.NoteOnEvent(tick=(time - lasttime) * lasting_on, channel=i, velocity=vel, pitch=one_pit))
                lasttime = time
                time = time + 1

        else:
            j = 1
            while j <= length - 1:
                if count_state >= len(state):
                    break
                one_state = state[count_state]
                if (length - 1 - j) < one_state:
                    break
                for k in range(one_state):
                    handled_pitch = handled_pitch + 1
                    if((handled_pitch % (960 / lasting)) == 0):
                        if handled_pitch < length:
                            track1 = get_track2(track1, pit[handled_pitch], music_combine)
                    one_pit = int(pit[j])
                    # one_vel = vel[j]
                    one_tick = tick[j]
                    if k == 0:
                        track.append(midi.NoteOnEvent(tick=(time - lasttime) * lasting_on, channel=i, velocity=vel, pitch=one_pit))
                        lasttime = time
                        time = time + 1
                    else:
                        track.append(midi.NoteOnEvent(tick=0, channel=i, velocity=vel, pitch=one_pit))
                    j = j + 1
                j = j - one_state
                for k in range(one_state):
                    handled_pitch = handled_pitch + 1
                    if((handled_pitch % (960 / lasting)) == 0):
                        if handled_pitch < length:
                            track1 = get_track2(track1, pit[handled_pitch], music_combine)
                    one_pit = int(pit[j])
                    # one_vel = vel[j]
                    one_tick = tick[j]
                    if k == 0:
                        track.append(midi.NoteOffEvent(tick=(time - lasttime) * lasting_off, channel=i, pitch=one_pit, velocity=vel))
                        lasttime = time
                        time = time + 1
                    else:
                        track.append(midi.NoteOffEvent(tick=0, channel=i, pitch=one_pit, velocity=vel))
                    j = j + 1
                count_state = count_state + 1
        eot = midi.EndOfTrackEvent(tick=2)
        track.append(eot)
        track1.append(eot)
    print pattern
    if not (os.path.exists("5_style{}".format(name1 + 1))):
        os.mkdir("5_style{}".format(name1 + 1))
    midi.write_midifile("5_style{}/style:{},songs:{}.mid".format(name1 + 1, name1 + 1, name2 + 1), pattern)
