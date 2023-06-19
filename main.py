import random as r
import numpy as np
import matplotlib.pyplot as plt

length = 10  # Размерность матрицы
length3 = 10  # Количество коэффициентов отделимости чисел
length4 = 10  # Количество дефформируемых матриц


# Создаем симметричную матрицу
def do_symmetric_matrix(l):
    a = [[1 for i in range(l)] for j in range(l)]
    l_ = int(l / 2)
    for i in range(l - l_):
        a[i][l - 1 - i] = r.randint(1, 10)
        a[l - 1 - i][i] = a[i][l - 1 - i]
    return a


a__ = np.array(do_symmetric_matrix(length))  # Получаем симметричную матрицу


# Создаем матрицу

def do_matrix(koef):
    mat_c = [[0 for i in range(length)] for j in range(length)]
    for i in range(length):
        mat_c[i][i] = (i * koef + 1)
        # mat_c[i][i] = 1 - i * koef
    mat_c = np.array(mat_c)  # Диагональная матрица
    # eig_value = np.array(result1[1])  # Собственные числа матрицы
    # eig_vector = np.array(result1[2])  # Матрица с точным значением собственных векторов
    # print('Заданые собственные числа: ')
    # print(eig_value)
    # print('Заданные собственные вектора: ')
    # print(eig_vector)
    a_ = np.array(np.linalg.inv(a__))  # Находим обратную матрицу
    matrix = a__.dot(mat_c).dot(a_)  # Получаем конечную симметричную матрицу
    return matrix


def deforming_matrix(a):
    massive_a = [[[0 for i in range(length)] for j in range(length)] for k in range(length4)]
    for k in range(length4):
        for i in range(length):
            for j in range(length):
                massive_a[k][i][j] = a[i][j] * (1 + 0.01 * k)
    return massive_a


# Реализуем метод плоских вращений Якоби


def Yakobi(a, e, n_max):
    current = [0 for j in range(length)]
    previous = [0 for k in range(length)]
    iteration = 0
    cur = 1
    eig_vector = np.eye(length)
    while cur > e and iteration < n_max:
        t = np.eye(length)
        e_m = 0
        for str in range(1, length):
            for col in range(str):
                if abs(a[str][col]) > abs(e_m):
                    e_m = a[str][col]
                    string = str
                    column = col
        p = 2 * e_m
        q = a[string][string] - a[column][column]
        if q == 0:
            c = (2 ** 0.5) / 2
            s = (2 ** 0.5) / 2
        else:
            d = (q ** 2 + p ** 2) ** 0.5
            rrr = abs(q) / (2 * d)
            c = (0.5 + rrr) ** 0.5
            s = ((0.5 - rrr) ** 0.5) * np.sign(p * q)
        t[string][string] = c
        t[column][column] = c
        t[string][column] = -s
        t[column][string] = s
        t_tranc = t.T
        a = t_tranc.dot(a).dot(t)
        if iteration == 0:
            eigenvectors_t = t
            for i in range(length):
                previous[i] = a[i][i]
            iteration += 1
            continue
        else:
            eigenvectors_t = eigenvectors_t.dot(t)
        for j in range(length):
            current[j] = abs(a[j][j] - previous[j])
            previous[j] = a[j][j]
        cur = max(current)
        if cur == 0:
            cur = cur2
            break
        cur2 = max(current)
        iteration += 1
    # Получаем окончательное собственное число
    eigenvalue = [0 for int2 in range(length)]
    for int in range(length):
        eigenvalue[int] = a[int][int]
    eigenvalue = eigenvalue[:: -1]
    # Получаем окончательный собственный вектор
    eigenvectors = [[0 for int5 in range(length)] for int3 in range(length)]
    sin_phi = [0 for i in range(length)]
    kin3 = 0
    for int3 in range(length):
        for kin2 in range(length):
            eigenvectors[kin3][kin2] = eigenvectors_t[int3][kin2]
        e_v = np.array(eigenvectors[kin3])
        phi = np.arccos((e_v.dot(eig_vector[int3])))
        sin_phi[int3] = abs(np.sin(phi))
        kin3 += 1
    rel_error = max(sin_phi)
    return rel_error, cur, iteration


# result2 = Yakobi(mat, 1e-16, 500)  # Записываем результат метода
# eig_value2 = result2[0]
# print('Получившиеся собственные числа:')
# print(eig_value2)
# eig_vector2 = result2[1]
# # print('Получившиеся собственные вектора:')
# # print(eig_vector2)
# relative_error = result2[2]
# # print('Синусы угла между векторами: ')
# # print(relative_error)
# # print('Относительная погрешность: ')
# # print(result2[3])
# # print('Количество итераций: ')
# # print(result2[4])


def chart1_2_3():
    coefficients = [0 for i in range(length3)]
    sinus = [0 for i in range(length3)]
    rel = [0 for i in range(length3)]
    iteration_max = [0 for i in range(length3)]
    for j in range(length3):
        coefficients[j] = 1 / (10 ** (j+1))
        mat_cf = do_matrix(coefficients[j])
        graph = Yakobi(mat_cf, 1e-16, 200)
        sinus[j] = graph[0]
        rel[j] = graph[1]
        iteration_max[j] = graph[2]
    return coefficients, sinus, rel, iteration_max


result2 = chart1_2_3()
coefficient = result2[0][::-1]
# print(coefficient)
sinus_vectors = result2[1]
# print(sinus_vectors)
relative_error = result2[2]
# print(relative_error)
iteration_max_ = result2[3]
# print(iteration_max_)


plt.figure(1)
plt.grid()
plt.xscale('log')
plt.xlabel('Коэффициент отделимости собственных чисел')
plt.ylabel('Синус угла между собственными векторами')
plt.title("График зависимости синуса угла между собственными векторами от отделимости собственных чисел")
plt.plot(coefficient, sinus_vectors)

plt.figure(2)
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Коэффициент отделимости собственных чисел')
plt.ylabel('Относительная погрешность')
plt.title("График зависимости относительной погрешности от отделимости собственных чисел")
plt.plot(coefficient, relative_error)

plt.figure(3)
plt.grid()
plt.xscale('log')
plt.xlabel('Коэффициент отделимости собственных чисел')
plt.ylabel('Максимальное количество итераций')
plt.title("График зависимости максимального количества итераций от отделимости собственных чисел")
plt.plot(coefficient, iteration_max_)


def chart4(p):
    rel = [0 for i in range(100)]
    iteration_max = [0 for i in range(100)]
    for j in range(100):
        iteration_max[j] = j + 101
        mat_cf = do_matrix(p)
        graph = Yakobi(mat_cf, 1e-16, iteration_max[j])
        rel[j] = graph[1]

    return rel, iteration_max


result3_1 = chart4(1e-1)
relative1 = result3_1[0]
number_iter1 = result3_1[1]

result3_2 = chart4(1e-5)
relative2 = result3_2[0]
number_iter2 = result3_2[1]

plt.figure(4)
plt.grid()
plt.yscale('log')
plt.xlabel('Номер итерации')
plt.ylabel('Относительная погрешность')
plt.title("График зависимости относительной погрешности от номера итерации")
plt.plot(number_iter1, relative1, label='Отделимость собственных чисел = 10')
plt.plot(number_iter2, relative2, label='Отделимость собственных чисел = 1.0001')
plt.legend()


def chart5(p):
    rel = [0 for i in range(10)]
    iteration_max = [0 for i in range(10)]
    for j in range(10):
        rel[j] = 1 / (10 ** (j + 3))
        mat_cf = do_matrix(p)
        graph = Yakobi(mat_cf, rel[j], 300)
        iteration_max[j] = graph[2]

    return rel, iteration_max


result4_1 = chart5(1e-1)
rel_err1 = result4_1[0]
iter_rel1 = result4_1[1]

result4_2 = chart5(1e-5)
rel_err2 = result4_2[0]
iter_rel2 = result4_2[1]

plt.figure(5)
plt.grid()
plt.xscale('log')
plt.xlabel('Задаваемая точность')
plt.ylabel('Максимальное количество итераций')
plt.title("График зависимости максимального количества итераций от задаваемой точности")
plt.plot(rel_err1, iter_rel1, label='Отделимость собственных чисел = 10')
plt.plot(rel_err2, iter_rel2, label='Отделимость собственных чисел = 1.0001')
plt.legend()


def chart6():
    rel = [0 for i in range(length4)]
    for j in range(length4):
        mat_cf = do_matrix(1e-1)
        mat_dm = np.array(deforming_matrix(mat_cf))
        graph = Yakobi(mat_dm[j], 1e-16, 300)
        rel[j] = abs(1 - (Yakobi(mat_dm[0], 1e-16, 300)[1] / graph[1]))
    return rel


error = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

plt.figure(6)
plt.grid()
plt.xlabel('Изменение исходной матрицы в процентах')
plt.ylabel('Относительная погрешность относительно исходной матрицы')
plt.title("График зависимости относительной погрешности от изменения исходной матрицы")
plt.plot(error, chart6(), label='Отделимость собственных чисел = 10')
plt.legend()
plt.show()
