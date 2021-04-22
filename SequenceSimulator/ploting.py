import numpy as np
import matplotlib.pyplot as plt


def ploting_function(structs):
    courses = []
    values = []
    i = int(0)
    name = ''
    for struct in structs:
        courses.append(structs[struct]['errors'])
        values.append(structs[struct]['accuracy'])
        name = structs[struct]['graph_name']
    courses = list(courses)
    values = list(values)
    fig = plt.figure(figsize=(10, 5))
    plt.bar(courses, values, color='maroon',
            width=0.4)

    plt.xlabel("Errors")
    plt.ylabel("Accuracy values")
    plt.title(name)
    plt.show()
