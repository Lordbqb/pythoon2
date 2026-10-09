# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:12:16 2026

@author: catko
"""
# %%
########## ПЕРЕМЕННЫЕ ##########
# Задание 1
stroka = "text"
integer = int(12)
floatnum = float(12)
logic = True
spisoc = [1, 2, 3]
kortech = (1, 2, 3)
mnozhestvo = {1, 2, 3}
slovar = {"a":1, "b":2, "c":3}
# вывод:
print(stroka, integer, floatnum,
logic, spisoc, kortech,
mnozhestvo, slovar, sep = '\n')

# %%
# Задание 2
list1 = [2, 3, 4, 5]
print(list1[0], list1[2])
list1.insert(0, 1)
list1.append(6)
print(list1, " len =", len(list1))
# %%
# Задание 3
dict1 = {"a":1, "b":2, "c":3}
print(dict1.keys(), dict1.values())
dict1["d"] = "new"
dict1["b"] = "update"
print(dict1)
# %%
# Задание 4
tuple1 = tuple(range(1,6))
print(tuple1[2], tuple1[4])
# кортеж индексирован, но не изменяем
# %%
# Задание 5
"""
Мне нравится тип boolean, потому что
в основе работы компьютеров лежит
смена 0 и 1 на транзисторах
"""
# %%
# Задание 6
(1+2 > 1) or (3*3 == 9) or ((7//2 != 3.5) and (not(2**2<5))) and (6%2 == 0) and (not(333/10 < 333))
# %%
########## ФУНКЦИИ ##########
# Задание 1
def F1(A, C, B=10):
    try:
        result = float((A + C) / B)
    except ZeroDivisionError:
        return "На ноль делить нельзя"
    return result
print(F1(2, 3))
print(F1(2, 3, 5))
print(F1(2, 3, 0))
# %%
# Задание 2
b = 100

def F2(a, c=None):
    if c is None:
        global b
        return a + b
    return a + c
print(F2(10))
print(F2(10, 5)) 
# %%
# Задание 3
def F3(a, b):
    if len(a) != len(b):
        raise ValueError("Строки разной длины")
    return {a[i]: b[i] for i in range(len(a))}
print(F3("abc", "123"))
# %%
# Задание 4
def F4(*args):
    if len(args) == 0:
        raise ValueError("Список аргументов пустой")
    return (sum(args) / len(args)) #сумма на кол-во элементов
print(F4(1, 2, 3, 4, 5, 6))
print(F4(8, 9, 10))
# %%
# Задание 5
def F51(a):
    return a**0.5
F51(9)

def F52(b):
    list52 = list()
    for e in b:
        list52.append(F51(e))
    return list52
F52([1, 4, 9, 16, 25])
# %%
########## ИСКЛЮЧЕНИЯ Try/except ##########
# Задание 1
def dev(a, b, c):
    try:
        a, b, c = map(float, (a, b, c))
    except (TypeError, ValueError):
        print("Аргументы не числа, операция не была проведена")
        return
    try:
        r = (a + b) / c
    except ZeroDivisionError:
        r = a + b
        comm = "c = 0, деление проведено не было:"
    else:
        comm = "Операция проведена успешно:"
    finally:
        print(comm, r)
        
dev("1", "1", "2")
dev("1", "1", "0")
dev("1", "abc", "2")
