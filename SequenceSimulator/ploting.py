import numpy as np
import matplotlib.pyplot as plt

# Function for plotting we needed for the presentation
def ploting_function(structs):
    courses = []
    values = []
    i = int(0)
    name = ''
    for struct in structs:
        courses.append(struct['errors'])
        values.append(struct['accuracy'])
        name = struct['graph_name']
    courses = list(courses)
    values = list(values)
    fig = plt.figure(figsize=(10, 5))
    plt.bar(courses, values, color='maroon',
            width=0.4)

    plt.xlabel("Errors")
    plt.ylabel("Accuracy values")
    plt.title(name)
    plt.show()
