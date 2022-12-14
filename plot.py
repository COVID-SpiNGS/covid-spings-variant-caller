import matplotlib.pyplot as plt


# Runtime per File
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y_python = [
    5.218637466430664,
    7.14712119102478,
    9.075253963470459,
    11.015209674835205,
    12.971837997436523,
    14.932650089263916,
    16.913827896118164,
    18.88459849357605,
    20.892345905303955,
    22.87979507446289
]
y_kotlin = [
    1.228,
    1.854,
    2.449,
    3.004,
    3.533,
    4.081,
    4.599,
    5.142,
    5.660,
    6.199
]

plt.ylabel('Runtime in seconds')
plt.xlabel('Number of files')
plt.plot(x, y_python, label='Python', color='#76B900')
plt.plot(x, y_kotlin, label='Kotlin', color='#0082D1')
plt.legend()
plt.show()

# Runtime per MB


x = [
    5950153 / 1024 / 1024,
    11871594 / 1024 / 1024,
    17796370 / 1024 / 1024,
    23712322 / 1024 / 1024,
    29621170 / 1024 / 1024,
    35534484 / 1024 / 1024,
    47366139 / 1024 / 1024,
    53274512 / 1024 / 1024,
    59194871 / 1024 / 1024
]
y_python = [
    5.295108795166016,
    7.579920530319214,
    10.147274255752563,
    12.721980571746826,
    16.16369652748108,
    19.27688694000244,
    26.451818227767944,
    30.67709183692932,
    34.758419036865234
]
y_kotlin = [
    1.053,
    1.141,
    1.379,
    1.589,
    1.723,
    1.875,
    2.249,
    2.540,
    2.601
]

plt.ylabel('Runtime in seconds')
plt.xlabel('Filesize in MB')
plt.plot(x, y_python, label='Python', color='#76B900')
plt.plot(x, y_kotlin, label='Kotlin', color='#0082D1')
plt.legend()
plt.show()


