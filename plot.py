import matplotlib.pyplot as plt


x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y_python = [
    5.485137224197388,
    4.183643102645874, 
    6.845597505569458,
    10.174527406692505,
    14.296257734298706,
    18.698223114013672,
    29.482377290725708,
    36.68915510177612,
    38.5833375453949,
    49.222522020339966
]
y_kotlin = [
    0.985,
    0.803,
    0.902,
    1.039,
    1.259,
    1.303,
    1.714,
    1.840,
    1.909,
    1.995
]

plt.ylabel('Runtime in seconds')
plt.xlabel('Number of files')
plt.plot(x, y_python, label="Python")
plt.plot(x, y_kotlin, label="Kotlin")
plt.legend()
plt.show()

