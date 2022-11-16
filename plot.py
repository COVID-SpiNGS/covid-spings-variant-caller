import matplotlib.pyplot as plt


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
plt.plot(x, y_python, label="Python")
plt.plot(x, y_kotlin, label="Kotlin")
plt.legend()
plt.show()


