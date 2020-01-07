# coding = utf-8
# encoding = utf-8
import tensorflow as tf
from pylab import *
from numpy import *
from math import sqrt
import random
import New_note_to_Midi as NM

mu = 0
sigma = 1
m_size = 16
d_epoches = 10
train_liters = 10


def conv2d(x, W):
    return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')


def max_pool(x):
    return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME')


def weight_variable(shape):
    initial = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(initial)


def bias_variable(shape):
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial)


def Conv_Net(input, w, h, z, out_dim, keep_prob=0.5):
    x_image = tf.reshape(input, [-1, w, h, z])
    W_conv1 = weight_variable([3, 3, z, 64])
    b_conv1 = bias_variable([64])
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool(h_conv1)
    _, m, n, _ = h_pool1.get_shape()
    h_pool1_new = tf.reshape(h_pool1, [1, int(m * n * 64)])
    W_fc1 = weight_variable([int(m * n * 64), out_dim])
    b_fc1 = bias_variable([out_dim])
    h_fc1 = tf.nn.sigmoid(tf.matmul(h_pool1_new, W_fc1) + b_fc1)
    out = tf.nn.dropout(h_fc1, keep_prob)
    return out, [W_conv1, b_conv1, W_fc1, b_fc1]


def Conv_Net_increase(input, out_dim):
    width = int(input.get_shape()[1])
    w1 = weight_variable([width, 8])
    b1 = bias_variable([8])
    out1 = tf.nn.relu(tf.matmul(input, w1) + b1)

    width = int(out1.get_shape()[1])
    w2 = weight_variable([width, out_dim])
    b2 = bias_variable([out_dim])
    out2 = tf.nn.relu(tf.matmul(out1, w2) + b2)
    return out2, [w1, b1, w2, b2]


def Discriminator(input, m, n):
    m = int(m)
    n = int(n)
    para = []
    conv1, para1_1 = Conv_Net(input, m, n, 1, 256)
    para.extend(para1_1)
    conv2, para1_2 = Conv_Net(conv1, 16, 16, 1, 64)
    para.extend(para1_2)
    conv3, para1_3 = Conv_Net(conv2, 8, 8, 1, 16)
    para.extend(para1_3)
    conv4, para1_4 = Conv_Net(conv3, 4, 4, 1, 1)
    para.extend(para1_4)
    conv5 = tf.reduce_mean(conv4)
    return conv5, para


def Generator(noise, z_size, prev_x, one_out_size):
    z_sqrt_size = int(sqrt(z_size))
    z_left_size = int(z_size / z_sqrt_size)
    out_sqrt_size = int(sqrt(one_out_size))
    out_left_size = int(one_out_size / out_sqrt_size)

    para = []
    pre_conv1, para1_1 = Conv_Net(prev_x, out_sqrt_size, out_left_size, 1, 256)
    para.extend(para1_1)
    pre_conv2, para1_2 = Conv_Net(pre_conv1, 16, 16, 1, 128)
    para.extend(para1_2)
    pre_conv3, para1_3 = Conv_Net(pre_conv2, 8, 16, 1, 64)
    para.extend(para1_3)
    pre_conv4, para1_4 = Conv_Net(pre_conv3, 8, 8, 1, 16)
    para.extend(para1_4)

    noise = tf.expand_dims(noise, 2)
    pre_conv4 = tf.expand_dims(pre_conv4, 2)
    input = tf.concat([noise, pre_conv4], axis=2)
    input_conv1, para2_1 = Conv_Net(input, 4, 4, 2, 64)
    para.extend(para2_1)
    input_conv1 = tf.expand_dims(input_conv1, 2)
    pre_conv3 = tf.expand_dims(pre_conv3, 2)
    input_conv1 = tf.concat([input_conv1, pre_conv3], axis=2)
    input_conv2, para2_2 = Conv_Net(input_conv1, 8, 8, 2, 128)
    para.extend(para2_2)
    input_conv2 = tf.expand_dims(input_conv2, 2)
    pre_conv2 = tf.expand_dims(pre_conv2, 2)
    input_conv2 = tf.concat([input_conv2, pre_conv2], axis=2)
    input_conv3, para2_3 = Conv_Net(input_conv2, 8, 16, 2, 256)
    para.extend(para2_3)
    input_conv3 = tf.expand_dims(input_conv3, 2)
    pre_conv1 = tf.expand_dims(pre_conv1, 2)
    input_conv3 = tf.concat([input_conv3, pre_conv1], axis=2)
    input_conv4, para2_4 = Conv_Net(input_conv3, 16, 16, 2, one_out_size)
    para.extend(para2_4)
    input_conv5 = tf.reshape(input_conv4, [1, one_out_size])
    return input_conv5, para


def optimizer(loss, var_list):
    batch = tf.Variable(0)
    learning_rate = tf.train.exponential_decay(
        0.01,
        batch,
        150,
        0.95,
        staircase=True
    )
    result = tf.train.AdamOptimizer(learning_rate).minimize(
        loss,
        global_step=batch,
        var_list=var_list
    )
    return result


def model(batch_size, one_out_size, times, position, pitch):
    length = np.shape(pitch[position])
    if length[0] >= batch_size * times:
        x_sqrt_size = int(sqrt(batch_size))
        out_sqrt_size = int(sqrt(one_out_size))
        # x_sqrt_size = sqrt(batch_size)
        with tf.variable_scope('G'):
            z_noi = tf.placeholder(tf.float32, shape=(1, m_size))
            pre = tf.placeholder(tf.float32, shape=(one_out_size, 1))
            G, theta_g = Generator(z_noi, m_size, pre, one_out_size)
            G = tf.multiply(G, 5.0)

        with tf.variable_scope('D') as scope:
            x_sam = tf.placeholder(tf.float32, shape=(batch_size, 1))
            D1, theta_d = Discriminator(x_sam, x_sqrt_size, int(batch_size / x_sqrt_size))
            D1 = tf.maximum(0.01, tf.minimum(0.99, D1))

            scope.reuse_variables()
            D2, theta_d = Discriminator(G, out_sqrt_size, int(one_out_size / out_sqrt_size))
            D2 = tf.maximum(0.01, tf.minimum(0.99, D2))

        loss_d = tf.reduce_mean(tf.log(D1) + tf.log(1 - D2))
        loss_g = tf.reduce_mean(tf.log(D2))

        opt_d = optimizer(-loss_d, theta_d)
        opt_g = optimizer(-loss_g, theta_g)

        sess = tf.InteractiveSession()
        tf.global_variables_initializer().run()

        histd, histg = np.zeros(train_liters), np.zeros(train_liters)

        z = np.linspace(mu - 5, mu + 5, m_size) + np.random.random(m_size) * 0.01
        pre_f = sess.run(G, {z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(pitch[position][0:batch_size], (batch_size, 1))})
        w = np.shape(pitch[position])
        max_epoches = int((w[0] / batch_size))
        epoch_lists = [x for x in range(max_epoches)]
        for i in range(train_liters):
            for j in range(d_epoches):
                choice = random.sample(epoch_lists, 1)[0] + 1
                x = pitch[position][(choice - 1) * batch_size:batch_size * choice]
                z = np.linspace(mu - 5, mu + 5, m_size) + np.random.random(m_size) * 0.01
                histd[i], _ = sess.run([loss_d, opt_d], feed_dict={x_sam: np.reshape(x, (batch_size, 1)), z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(pre_f, (one_out_size, 1))})
                pre_f = sess.run(G, {z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(pre_f, (one_out_size, 1))})
            z = np.linspace(mu - 5, mu + 5, m_size) + np.random.random(m_size) * 0.01
            histg[i], _ = sess.run([loss_g, opt_g], {z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(pre_f, (one_out_size, 1))})

        pit_all = []
        pit_max = max(pitch[position])
        pit_min = min(pitch[position])
        for i in range(times):
            z = np.linspace(mu - 5, mu + 5, m_size)
            if i == 0:
                g_new = sess.run(G, {z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(pitch[position][0:one_out_size], (one_out_size, 1))})
                each_out = g_new
            else:
                g_new = sess.run(G, {z_noi: np.reshape(z, (1, m_size)), pre: np.reshape(each_out, (one_out_size, 1))})
                each_out = g_new
            one_git = each_out[0].astype(int)
            # print(each_out[0])
            tran = np.fft.fft(one_git)
            filter = np.concatenate((np.ones(200), np.zeros(200)))
            one_git = np.fft.ifft(tran * filter)
            new_ti = (one_git - min(one_git)) / (max(one_git) - min(one_git)) * (pit_max - pit_min) + pit_min
            new_ti = new_ti.astype(int)
            new_ti = new_ti.tolist()
            pit_all.extend(new_ti)
        return pit_all
    else:
        return 0
