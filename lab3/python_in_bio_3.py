# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:17:36 2026

@author: Kokhanova E. A.
"""
# %%
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
# %%
# ГЛОБАЛЬНОЕ ВЫРАВНИВАНИЕ
seq1 = Seq("AAGCACCGGACCAGGGGGAAGGCGAGCAGTGCCAATCTACAGCGAAGAAAGTCTCGTTTGGTAAAAGCGAGAGGGGAAAGCCTGAGCATGCAGAGTGTGCAGAGCACGAGCTTTTGTCTCCGAAAGCAGTGCCTTTGCCTGACCTTCCTGCTTCTCCATCTCCTGGGACAGGTAAGTGGCACACCCTTAAGATGCCCCCAAGTTACTTTGCCCGCCTTGGTGGCCCCCATTTGGTCACCGGGCTCACTGCGTCTTCTGTCCCAGCTGAGTGGTTTCTCCTTGTCTCGCCTGCCTTCAGGTCGCTGCGACTCAGCGCTGCCCTCCCCAGTGCCCGGGCCGGTGCCCTGCGACGCCGCCGACCTGCGCCCCCGGGGTGCGCGCGGTGCTGGACGGCTGCTCATGCTGTCTGGTGTGTGCCCG") 
seq2 = Seq("AAGCACCGGACCAGGGGGAAGGCGAGCAGTGCCAATCTACAGCGAAGAAAGTCTCGTTTGGTAAAAGCGAGAGGGGAAAGCCTGAGCATGCAGAGTGTGCAGAGCACGAGCTTTTGTCTCCGAAAGCAGTGCCTTTGCCTGACCTTCCTGCTTCTCCATCTCCTGGGTAAGTGACAGGGCACACCCTTAAGATGCCCCCAAGTTACTTTGCCCGCCTTGGTGGCCCCCATTTGGTCACCGGGTGCGTCTCTCACTCTGTCCCAGCTGAGTGGTTTCTCCTTGTCTCGCCTGCCTTCAGGTCGCTGCGACTCAGCGCTGCCCTCCGCCCGGCCAGTGCCGGTGCCCTGCGACGCCGCCGACCTGCGCCCCCGGGGTGCGCGCGGTGCTGGACGGCTGCTCATGCTGTCTGGTGTGTGCCCG")

# Создаем объект выравнивателя
aligner1 = PairwiseAligner()

# Устанавливаем режим (глобальное выравнивание)
aligner1.mode = "global" 

# Задаем параметры (вместо знаков вопроса из задания)
# балл за совпадение
aligner1.match_score = 2
# штраф за несовпадение
aligner1.mismatch_score = -1
# штрафы за пропуски/разрывы
aligner1.open_gap_score = -5   #— штраф за открытие гэпа
aligner1.extend_gap_score = -1 #— штраф за продолжение гэпа

# Выполняем выравнивание
res_align = aligner1.align(seq1, seq2)
len(res_align)
print(len(res_align), res_align[0], sep="\n")
# %%
# СОХРАНЕНИЕ В FASTA ФАЙЛ

# print(res_align[0][0], res_align[0][1], sep="\n")
# print(len(res_align[0][0]), len(res_align[0][1]), sep=" ; ")
# type(res_align[0][0])

# Создаем запись для первой последовательности
record1 = SeqRecord(Seq(res_align[0][0]), id="Origin1", name= "seq1 (aligned)", description="")

# Создаем запись для второй последовательности
record2 = SeqRecord(Seq(res_align[0][1]), id="Origin2", name= "seq2 (aligned)", description="")

print(record1, record2, sep="\n\n")
records = [record1, record2]
SeqIO.write(records, "result.fa", "fasta")
# %%
# ЛОКАЛЬНОЕ ВЫРАВНИВАНИЕ
seq3 = "CCATCTCCTGGGACA"

aligner2 = PairwiseAligner()

aligner2.mode = "local" 

aligner2.match_score = 2
aligner2.mismatch_score = -1
aligner2.open_gap_score = -5
aligner2.extend_gap_score = -1

loc_align = aligner2.align(seq1, seq3)
print(len(loc_align), loc_align[0], sep="\n")
# %%

# %%
# ПРОСТОЕ ХЭШИРОВАНИЕ
def my_hash(text):
    return sum(ord(char) for char in text)
print([my_hash(word) for word in ['cat', 'dog', 'god', 'thing', 'night']])
## для слов с одинаковыми буквами хэш одинаковый##
# %%
# ХЭШ-ТАБЛИЦА СЛОВ
def get_word_frequencies(text):
    # приводим к нижнему регистру и меняем знаки препинания на пробел
    clean_text = text.lower().replace(",", "").replace(".", "").replace("’", "")
    words = clean_text.split()
    
    freq_dict = {}
    for w in words:
        if w in freq_dict:
            freq_dict[w] += 1
        else:
            freq_dict[w] = 1
    return freq_dict

str1 = 'I felt happy because I saw the others were happy and because I knew I should feel happy, but I wasn’t really happy.'

word_counts = get_word_frequencies(str1)
print("Словарь частот:", word_counts)

most_frequent = max(word_counts, key=word_counts.get)
print("Самое частое слово:", most_frequent)

repeated_words = [word for word, count in word_counts.items() if count > 1]
print("Слова, встречающиеся > 1 раза:", repeated_words)
# %%
# ХЭШ-ТАБЛИЦА ПОДСТРОК
dna = "ATCAGGCATCAGGCACATGATCGTACTACATGAATCAGTGAATCGATCAGTACACGT"

def build_hash_table(seq, k=4):
    hash_table = {}
    
    # Проходим по всей строке. 
    # len(seq) - k + 1 — это количество возможных окон.
    for i in range(len(seq) - k + 1):
        fragment = seq[i:i+k] # Вырезаем окно длиной k
        if fragment in hash_table:
            hash_table[fragment].append(i) # если фрагмент уже есть, добавляем позицию
        else:
            hash_table[fragment] = [i]  # если нет, создаем новый список
            
    return hash_table


dna_table = build_hash_table(dna, k=4)

print(f"Позиции подстроки 'ATCA': {dna_table.get('ATCA')}")