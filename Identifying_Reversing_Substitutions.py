# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 16:02:27 2026

@author: catko
"""

# у нас дано дерево, узлы которого — вид и его последовательность

# %%
# создадим объект, в атрибутах которого записано: 
# # положение узла (предок и потомки),
# название вида,
# последовательность вида:

class Node:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parent = None
        self.seq = None
    
    def __str__(self): #для красивого вывода
        return f"{self.__class__.__name__}({self.name})"
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.name})"
# %%
# начнём считывать датасет с того, что разделим структуру дерева и FASTA

with open("Sample Dataset.txt", "r") as f:
    lines = f.readlines()

tree_string = lines[0].strip().rstrip(";") # первая строка со структурой дерева

fasta_strings = [] # список с очищенными строками, в которых записаны FASTA
for line in lines[1:]:
    fasta_strings.append(line.strip())

# %%
# обработаем FASTA
sequences = {}
current_name = None
for line in fasta_strings:
    if line[0] == ">":
        current_name = line[1:]
    else:
        sequences[current_name] = line
# теперь sequences это словарь, 
# где кадому виду стоит в соответствие последовательность
# %%
def parse_newick(tree_string):
    # ============================================================
    # ПОДГОТОВКА
    # ============================================================
    
    s = tree_string.strip()
    n = len(s)
    i = 0
    stack = []
    root = None
    
    # ============================================================
    # ВСПОМОГАТЕЛЬНАЯ ФУНКЦИЯ: читаем имя вида
    # ============================================================

    def read_name(i):
        """
        Начиная с позиции i, читает символы, пока не наткнётся 
        на скобку, запятую или конец строки.
        Возвращает: (имя, новая_позиция_i)
        """
        start = i
        while i < n and s[i] not in "(,)":
            i += 1
        return s[start:i], i
    
    # ============================================================
    # ГЛАВНЫЙ ЦИКЛ: идём по строке слева направо
    # ============================================================
    
    while i < n:
        c = s[i]
        
        # --------------------------------------------------------
        # СЛУЧАЙ 1: открывающая скобка "("
        # --------------------------------------------------------
        if c == "(":
            # Открылась новая группа (внутренний узел).
            node = Node(name=None)
            
            # Если стек не пуст — значит, у этого узла есть родитель
            if stack:
                stack[-1].children.append(node)   # родитель запоминает потомка
                node.parent = stack[-1]           # потомок запоминает родителя
            
            # Кладём новый узел в стек
            stack.append(node)
            i += 1      # сдвигаем курсор за скобку
        
        # --------------------------------------------------------
        # СЛУЧАЙ 2: запятая ","
        # --------------------------------------------------------
        elif c == ",":
            i += 1

        # --------------------------------------------------------
        # СЛУЧАЙ 3: закрывающая скобка ")"
        # --------------------------------------------------------
        elif c == ")":
            # Закрылась группа, узел на вершине стека полностью сформирован
            closed = stack.pop()    # достаём его из стека
            
            # Сразу после ")" идёт имя этого узла — читаем его
            i += 1
            name, i = read_name(i)
            
            if name:
                closed.name = name          # даём узлу имя
            
            # Если стек опустел после pop — значит, мы только что закрыли
            # самый верхний узел, у которого нет родителя
            if not stack:
                root = closed
        
        # --------------------------------------------------------
        # СЛУЧАЙ 4: буква (начало имени)
        # --------------------------------------------------------
        else:
            # Мы встретили имя, которое идёт не после ")" —
            # это лист (конечный вид)
            name, i = read_name(i)
            
            leaf = Node(name=name)   # лист сразу создаём с именем
            
            if stack:
                # привязываем лист к тому, кто сейчас на вершине стека
                stack[-1].children.append(leaf)
                leaf.parent = stack[-1]
            else:
                # стек пуст — значит, это единственный узел в дереве
                root = leaf

    # ============================================================
    # ВОЗВРАЩАЕМ КОРЕНЬ
    # ============================================================
    return root        
# %%
# Проверим работу функции
def print_tree(node, level=0):
    print("  " * level + f"- {node.name}")
    for child in node.children:
        print_tree(child, level + 1)

root = parse_newick("((dog,cat)pet,fox)animal")
print_tree(root)

print("----------")

root = parse_newick(tree_string)
print_tree(root)
# %%
# Заполним seq у видов:
def fill_sequences(node, sequences):
    
    if node is None:
        return
    
    if node.name in sequences:
        node.seq = sequences[node.name]
    else:
        node.seq = None
    
    # рекурсивно обрабатываем всех потомков
    for child in node.children:
        fill_sequences(child, sequences)
# %%
# Проверим работу функции
def print_tree_with_seq(node, level=0):
    print("  " * level + f"- {node.name}: {node.seq}")
    for child in node.children:
        print_tree_with_seq(child, level + 1)

root = parse_newick(tree_string)
fill_sequences(root, sequences)
print_tree_with_seq(root)
# %%
# ПЕРЕХОДИМ К ПОИСКУ РЕВЕРСИЙ
# %%
# Создадим список всех видов, описанных в дереве, чтобы
# было удобно перебирать виды при просмотре их
# родителей и детей
def get_all_nodes(root):
    if root is None:
        return []
    
    result = [root]
    
    for child in root.children: # применяем рекурсию
        result.extend(get_all_nodes(child)) 
        
    return result

root = parse_newick(tree_string)
get_all_nodes(root)
# %%
# найдём мутации (по одной, указанной позиции)
def find_mutation_edges(all_nodes, i):
    edges = []
    for species in all_nodes:
        p = species.parent
        if p is None:
            continue                      # корень — у него нет родителя
        if p.seq[i] != species.seq[i]:    # символы различаются → мутация
            edges.append((p, species))
    return edges
# %%
# выведем мутации
root = parse_newick(tree_string)
fill_sequences(root, sequences)
def print_mutation(root, i):
    edges = find_mutation_edges(get_all_nodes(root), i)
    for ori, mut in edges:
        print(f"{ori.name} -> {mut.name}  ({ori.seq[i]} -> {mut.seq[i]})")

print_mutation(root, 0)
# %%
# поиск реверсии по одной мутации:
def find_reversions_from_edge(mutation_child, position, original_char, new_char):
    """
    mutation_child — узел-потомок на ребре-мутации (у него символ = new_char)
    position       — индекс позиции в последовательности
    original_char  — символ, который был у родителя (до мутации)
    new_char       — символ, который стал у потомка (после мутации)
    """
    reverting_nodes = []
    stack = [mutation_child]
    
    while stack:
        current_node = stack.pop()
        
        for child in current_node.children:
            child_char = child.seq[position]
            
            if child_char == original_char:
                reverting_nodes.append(child)
            
            elif child_char == new_char:
                stack.append(child)
            
            # третий случай (child_char == что-то другое) игнорируем
    
    return reverting_nodes
# %%
# Проверим работу функции
i = 0
edges = find_mutation_edges(get_all_nodes(root), i)

for mutation_parent, mutation_child in edges:
    original_char = mutation_parent.seq[i]
    new_char = mutation_child.seq[i]
    
    print(f"Ребро-мутация: {mutation_parent.name} -> {mutation_child.name}")
    print(f"позиция мутации: {i+1}")
    print(f"замена: {original_char} -> {new_char}")
    
    reverting_nodes = find_reversions_from_edge(
        mutation_child, i, original_char, new_char)
    
    if not reverting_nodes:
        print("  реверсий нет")
    else:
        for node in reverting_nodes:
            print(f"  реверсия: {mutation_child.name} {node.name} {i+1} "
                  f"{original_char}->{new_char}->{original_char}")
    print("")
# %%
# выводим все мутации
root = parse_newick(tree_string)
fill_sequences(root, sequences)
all_nodes = get_all_nodes(root)

seq_length = len(root.seq)   # длина последовательности одинакова у всех узлов
results = []

for i in range(seq_length):
    edges = find_mutation_edges(all_nodes, i)
    
    for mutation_parent, mutation_child in edges:
        
        original_char = mutation_parent.seq[i]
        new_char = mutation_child.seq[i]
        
        reverting_nodes = find_reversions_from_edge(
            mutation_child, i, original_char, new_char
        )
        
        for node in reverting_nodes:
            results.append(
                f"{mutation_child.name} {node.name} {i+1} "
                f"{original_char}->{new_char}->{original_char}"
            )

# печатаем все найденные реверсии
for line in results:
    print(line)
# %%
