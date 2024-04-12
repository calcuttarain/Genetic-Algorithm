from math import log2, ceil
import random
import matplotlib.pyplot as plt

#citirea fisierului de intrare
def read_file(file_path):
    global no_chromosomes, domain, coefficients, precision, crossover_probability, mutation_probability, no_generations
    file = open(file_path, 'r')
    no_chromosomes = int(file.readline())
    domain = [int(x) for x in file.readline().split()]
    coefficients = [int(x) for x in file.readline().split()]
    precision = int(file.readline())
    crossover_probability = float(file.readline())
    mutation_probability = float(file.readline())
    no_generations = int(file.readline())
    discretize()
    

    
def write_file(initial_generation, mean, mx, x_mx):
    file = open("output.txt", "w")
    encoded = encode(initial_generation)
    fitness = list(map(fitness_function, initial_generation))
    selection_chance = chances(initial_generation)
    file.write("Populatia initiala:\n")
    for i in range(len(initial_generation)):
        file.write(f"Cromozomul {i+1}: {encoded[i]}, x = {initial_generation[i]}, f(x) = {fitness[i]}\n")
    file.write("\nProbabilitati de selectie pentru populatia initiala:\n")
    for i in range(len(selection_chance)):
        file.write(f"Cromozomul {i+1}: {selection_chance[i]}\n")
    intervals_selection_chance = selection_chances_intervals(selection_chance)
    file.write("\nIntervale probabilitati selectie:\n")
    for i in intervals_selection_chance:
        file.write(str(i) + " ")
    file.write("\n\n")
    p1 = []
    for i in range(len(selection_chance)):
        u = random.random()
        index = find_interval(intervals_selection_chance, u)
        p1.append(initial_generation[index])
        file.write(f"u = {u} => selectam cromozomul numarul {index + 1}\n")
    file.write("\nDupa selectie:\n")
    encoded = encode(p1)
    fitness = list(map(fitness_function, p1))
    for i in range(len(p1)):
        file.write(f"Cromozomul {i+1}: {encoded[i]}, x = {p1[i]}, f(x) = {fitness[i]}\n")
        
    file.write("\nCine participa la cross-over:\n")
    indici = []
    for i in range(len(p1)):
        u = random.random()
        if u < crossover_probability:
            indici.append(i)
            file.write(f"u = {u} < {crossover_probability} participa\n")
        else:
            file.write(f"u = {u}\n")
    if len(indici) % 2 == 1:
        indici.pop()
    i = 0
    while i != len(indici):
        file.write(f"Cross-over intre cromozomul {indici[i]+1} si cromozomul {indici[i+1]+1}:\n")
        taietura = random.randint(0, chromosome_length)
        c1 = encoded[indici[i]]
        c2 = encoded[indici[i + 1]]
        prel_c1 = c2[:taietura] + c1[taietura:]
        prel_c2 = c1[:taietura] + c2[taietura:]
        encoded[indici[i]] = prel_c1
        encoded[indici[i + 1]] = prel_c2
        file.write(f"{prel_c1}, {prel_c2}, taietura: {taietura}\n")
        i += 2
    p1 = decode(encoded)
    fitness = list(map(fitness_function, p1))
    file.write("\nDupa cross-over:\n")
    for i in range(len(p1)):
        file.write(f"Cromozomul {i+1}: {encoded[i]}, x = {p1[i]}, f(x) = {fitness[i]}\n")
    
    file.write("\nCromozomii care vor participa la mutatie: ")
    for i in range(len(encoded)):
        u = random.random()
        print(u)
        if u < mutation_probability:
            print("da")
            taietura = random.randint(0, chromosome_length)
            file.write(f"{i+1} ")
            mutagen = ""
            for bit in encoded[i][:taietura]:
                if bit == "0":
                    mutagen += "1"
                else:
                    mutagen += "0"
            mutagen += encoded[i][taietura:]
            encoded[i] = mutagen
    
    p1 = decode(encoded)
    fitness = list(map(fitness_function, p1))
    file.write("\nDupa mutatie:\n")
    for i in range(len(p1)):
        file.write(f"Cromozomul {i+1}: {encoded[i]}, x = {p1[i]}, f(x) = {fitness[i]}\n")
    file.write("\nPentru restul generatiilor, media fitness vs maximul fitness:\n")
    for i in range(len(mx)):
        file.write(f"{mean[i]} {mx[i]}\n")
    evolutie(mx, x_mx)
    
def maxim_functie_grad2(a, b, c, interval):
    f_a = a * interval[0] ** 2 + b * interval[0] + c
    f_b = a * interval[1] ** 2 + b * interval[1] + c

    x_crit = -b / (2 * a)
    if interval[0] <= x_crit <= interval[1]:
        f_crit = a * x_crit ** 2 + b * x_crit + c
    else:
        f_crit = float('-inf')
    maxim = max(f_a, f_b, f_crit)
    if maxim == f_a:
        x_maxim = interval[0]
    elif maxim == f_b:
        x_maxim = interval[1]
    else:
        x_maxim = x_crit

    return x_maxim, maxim
    
def evolutie(mx, x_mx):
    plt.figure(figsize=(8, 6))

    plt.plot(x_mx, mx, marker='o', color='blue', linestyle='-', label='f(x)')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    x, y = maxim_functie_grad2(coefficients[0], coefficients[1], coefficients[2], domain)
    plt.scatter([x], [y], color='red', label='maximul lui f(x) (aproximat)')
    plt.title('Evolutia algoritmului genetic')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    
#calcularea numarului de biti pentru fiecare cromozom, a pasului de discretizare si a intervalelor
def discretize():
    global chromosome_length, discretization_step, intervals
    chromosome_length = ceil(log2((domain[1] - domain[0]) * (10 ** precision)))
    discretization_step = (domain[1] - domain[0])/(2 ** chromosome_length)
    intervals = [[domain[0] + i * discretization_step, domain[0] + (i + 1) * discretization_step] for i in range (2 ** chromosome_length - 1)]

#cautarea binara pentru incadrarea unui cromozom in intervalul din care face parte din lista de intervale
def binary_search_intervals(x):
    left = 0
    right = len(intervals) - 1
    while left <= right:
        mid = (left + right) // 2
        interval = intervals[mid]
        if interval[0] <= x < interval[1]:
            return mid
        elif x < interval[0]:
            right = mid - 1
        else:
            left = mid + 1
    
#codificarea cromozomilor
def encode(numbers):
    encoded_numbers = []
    for number in numbers:
        index = binary_search_intervals(number)
        binary_index = format(index, 'b')
        encoded_number = '0' * (chromosome_length - len(binary_index)) + binary_index
        encoded_numbers.append(encoded_number)
    return encoded_numbers

#decodificarea cromozomilor(capatul inferior al intervalului din care face parte)
def decode(numbers):
    decoded_numbers = []
    for number in numbers:
        index = int(number, 2)
        decoded_numbers.append(intervals[index][0])
    return decoded_numbers

def fitness_function(x):
    return coefficients[0] * x**2 + coefficients[1] * x + coefficients[2]

def selection_chances_intervals(selection_chance):
    intervals_selection = []
    partial = 0
    for i in range(len(selection_chance)):
        intervals_selection.append(partial)
        partial += selection_chance[i]
    intervals_selection.append(round(partial, 2))
    return intervals_selection

def chances(generation):
    p = list(map(fitness_function, generation))
    total = sum(p)
    p = list(map(lambda x: x/total, p))
    return p

def find_interval(numbers, u):
    n = len(numbers)
    left = 0
    right = n - 1

    while left < right:
        mid = (left + right) // 2
        if numbers[mid] <= u:
            left = mid + 1
        else:
            right = mid
            
    if numbers[left] <= u:
        return left
    else:
        return left - 1

def selection(generation):
    selection_chance = chances(generation)
    elitist = generation[selection_chance.index(max(selection_chance))]
    selection_chances_interval = selection_chances_intervals(selection_chance)
    p1 = []
    for i in range(len(selection_chance)):
        u = random.random()
        index = find_interval(selection_chances_interval, u)
        p1.append(generation[index])
    if elitist not in p1:
        p1.pop()
        p1.append(elitist)
    return p1

def cross_over(sir):
    indici = []
    for i in range(len(sir)):
        u = random.random()
        if u < crossover_probability:
            indici.append(i)
    if len(indici) % 2 == 1:
        indici.pop()
        
    i = 0
    while i != len(indici):
        taietura = random.randint(0, chromosome_length)
        c1 = sir[indici[i]]
        c2 = sir[indici[i + 1]]
        prel_c1 = c2[:taietura] + c1[taietura:]
        prel_c2 = c1[:taietura] + c2[taietura:]
        sir[indici[i]] = prel_c1
        sir[indici[i + 1]] = prel_c2
        i += 2
    return sir
        
def mutation(sir):
    for i in range(len(sir)):
        u = random.random()
        if u < mutation_probability:
            taietura = random.randint(0, chromosome_length)
            mutagen = ""
            for bit in sir[i][:taietura]:
                if bit == "0":
                    mutagen += "1"
                else:
                    mutagen += "0"
            mutagen += sir[i][taietura:]
            sir[i] = mutagen
    return sir
        
def genetic_algorithm(file_path):
    read_file(file_path)
    generation = [round(random.uniform(domain[0], domain[1]), precision) for _ in range(no_chromosomes)]
    initial_generation = generation
    n = no_generations
    mean = []
    mx = []
    x_mx = []
    while n:
        p1 = selection(generation)
        p1 = encode(p1)
        p2 = cross_over(p1)
        generation = mutation(p2)
        generation = decode(generation)
        
        p = list(map(fitness_function, generation))
        mx.append(max(p))
        x_mx.append(generation[p.index(max(p))])
        mean.append(sum(p)/len(p))
        
        n -= 1

    write_file(initial_generation, mean, mx, x_mx)
    
genetic_algorithm("input.txt")

