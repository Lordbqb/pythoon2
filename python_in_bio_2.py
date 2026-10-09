# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 12:10:53 2026

@author: Kokhanova E A

"""
# %%
########## КЛАССЫ ##########
# Задание 1.1
from abc import ABC, abstractmethod

### абстрактный класс "девайс" ###
class Device(ABC):
    # список разрешённых типов девайсов
    _allowed_subclasses = frozenset({"Lamp", "Thermostat", "Socket"})

    def __init_subclass__(device_class, **kwargs):
        super().__init_subclass__(**kwargs)
        if device_class.__name__ not in Device._allowed_subclasses:
            raise TypeError(
                f"Класс {device_class.__name__} запрещён: девайсы ограничены списком "
                f"{sorted(Device._allowed_subclasses)}"
            )

    def __init__(self, name):
        self.name = name
        self._is_on = False

    @abstractmethod
    def turn_on(self):
        ...

    @abstractmethod
    def turn_off(self):
        ...
    
    @property
    def is_on(self):
        return self._is_on

    def status(self):
        return f"{self.name}: {'ON' if self._is_on else 'OFF'}"

    def __str__(self): #для красивого вывода
        return f"{self.__class__.__name__}({self.name})"
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.name})"


### конкретные типы "девайсов" ###
class Lamp(Device):
    def turn_on(self):
        self._is_on = True
        return f"Лампа {self.name} светит"

    def turn_off(self):
        self._is_on = False
        return f"Лампа {self.name} погасла"

class Thermostat(Device):
    # переопределение __init__, чтобы добавить параметр температуры
    def __init__(self, name, temp=22):
        super().__init__(name)
        self._temp = temp

    def turn_on(self):
        self._is_on = True
        return f"Термостат {self.name} нагревает до {self._temp}°C"

    def turn_off(self):
        self._is_on = False
        return f"Термостат {self.name} выключен"

    def set_temperature(self, value):
        if not 10 <= value <= 40:
            raise ValueError("Температура должна быть от 10 до 40")
        self._temp = value
        return f"Температура установлена на {value}°C"

class Socket(Device):
    def turn_on(self):
        self._is_on = True
        return f"Розетка {self.name} под напряжением"

    def turn_off(self):
        self._is_on = False
        return f"Розетка {self.name} обесточена"


### класс "умного дома", содержащий в себе "девайсы" ###
class SmartHouse:
    def __init__(self, address):
        self.address = address
        self._temperature = 22
        self.__devices = [] #приватный, чтобы нельзя было случайно изменить

    def add_device(self, device):
        if not isinstance(device, Device):
            raise TypeError("Можно добавлять только девайсы")
        self.__devices.append(device)

    def get_devices(self):
        return list(self.__devices)

    def turn_on_all(self):
        return [d.turn_on() for d in self.__devices]

    def turn_off_all(self):
        return [d.turn_off() for d in self.__devices]

    def set_temperature(self, value):
        if not 10 <= value <= 40:
            raise ValueError("Температура должна быть от 10 до 40")
        self._temperature = value

    def get_temperature(self):
        return self._temperature

    def status(self):
        return {
            "address": self.address,
            "temperature": self._temperature,
            "devices": [d.status() for d in self.__devices],
            "devices_on": [d.name for d in self.__devices if d.is_on]
        }

    def __str__(self):
        return f"SmartHouse({self.address}, устройств: {len(self.__devices)})"

# %%
### проверки задания 1###
Device("Абстрактное устройство")

class Wunderwaffe(Device):
    def turn_on(self): return "Бдыщ"
    def turn_off(self): return "..."
    
house = SmartHouse("ул. Ленина, 10")
house.add_device(Lamp("Настольная"))
house.add_device(Thermostat("Кухонный", 24))
house.add_device(Socket("У окна"))

print(house.turn_on_all())
print(house.turn_off_all())

house.add_device("лампа")
house.add_device(42)

house.__devices
house.get_devices()

house.set_temperature(50)
house.set_temperature(25)
house.get_devices()[0].turn_on()
print(house.status())


# %%
########## МАТРИЧНЫЕ ПРОЦЕССЫ ##########
#import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.MeltingTemp import Tm_NN

# %%
# Задание 2.1
raw_seq =  '''        1 ctgactcaaa cgtctcactg atcgttacac aacaccagtg taatttcgca gtgcttatta
       61 acaagctaaa gcgctagcat ttactggaag cacgacttcc tctcgccaac gtggatataa
      121 agctgtatta atcatgggtg ttcacgactt gtggtccatc gtggagccag tccgagaggc
      181 ggtgccgttg tacagtttga gcggaaagac tctggcggta gacctgagct tatgggtctg
      241 tgaggctcag cAcgtccagg ccatgatggg cagagttatg aagccacatc tgaggaattt
      301 atttttccga gtatcatcgc tcacactgat GGGGtgaag ctcgtctttg tcatggaggg
      361 agaagctccc caacttaaat cagaaacaat gagcaaaagg acagaaatga gattcggagg
      421 atccagaaaa gCCcatttct caaagtcatc cacaactacc agcagagggc gctttaaagc
      481 tgtacttcga gGGGtgtgcag atatgttgaa ctacctgggt gtgccctggg tgacagcagc
      541 tggggaagca gaggccatgt gtgcttacct ggactcagca gggctggttg atggctgcat
      601 caccaacgat ggagatgcat tcttatatgg ggctaggacg gtctacagga acttcaatat
      661 gaacagtaaa gaccctcagg tcgattgcta ccagacatcc cgagtgcaga ccgagctaca
      721 cctgtcaaga gagaatcttg ttggcttggc cattctcctc ggttgtgatt aagcagatgt
      781 gttgtataaa gtaat'''
      
def intriplet(raw_seq, raw_output = True):
    seq = ''
    for ch in raw_seq:
        ch = ch.upper()
        if ch in 'ATGC':
            seq += ch
    
    if raw_output:
        triplets = "".join([seq[i:i+3] for i in range(0, len(seq), 3)])
    else:
        triplets = " ".join([seq[i:i+3] for i in range(0, len(seq), 3)])

    return triplets

intriplet(raw_seq, False)
# %%
# Задание 2.2 (через биопайтон)
intriplet(raw_seq)
dna = Seq(intriplet(raw_seq))

print(f"Длина цепи: {len(dna)} нуклеотидов")
print(f"GC-состав: {(gc_fraction(dna) * 100):.2f}%")
print(f"Молекулярная масса: {molecular_weight(dna, seq_type='DNA'):.2f} г/моль")
print(f"Температура плавления: {Tm_NN(dna) :.2f} °C")
# %%
# Задание 2.2 (через чистый пайтон)
clean_seq = intriplet(raw_seq)

length = len(clean_seq)

gc_count = clean_seq.count('G') + clean_seq.count('C')
gc_content = gc_count / length * 100

weights = {'A': 313.21, 'T': 304.20, 'G': 329.21, 'C': 289.18}
mass = sum(weights[base] for base in clean_seq) + 18.015 #+масса воды

#Tm = 64.9 + 41 * (кол-во G и C - 16.4) / длину
tm = 64.9 + 41 * (gc_count - 16.4) / length

print(f"Длина цепи: {length} нуклеотидов")
print(f"GC-состав: {gc_content:.2f}%")
print(f"Молекулярная масса: {mass:.2f} г/моль")
print(f"Температура плавления (прибл.): {tm:.2f} °C")
# %%
# Задание 2.3 (через биопайтон)
intriplet(raw_seq)
dna = Seq(intriplet(raw_seq))
mrna = dna.transcribe()
print(mrna)
# %%
# Задание 2.3 (через чистый пайтон)
mrna = intriplet(raw_seq).replace('T', 'U')
print(mrna)
# %%
# Задание 2.4
mrna = intriplet(raw_seq).replace('T', 'U')
def mrna_to_protein(mrna):
    try:
        start = str(mrna).index('AUG')
    except ValueError:
        raise ValueError("В последовательности нет старт-кодона AUG")
    
    coding_seq = Seq(str(mrna)[start:])
    protein = coding_seq.translate(to_stop=True) #до первого стоп-кодона
    
    return list(protein)

amino_acids = mrna_to_protein(str(mrna))
print(amino_acids)
# %%
# Задание 2.5
import random

class Protein:
    """Класс для работы с белковыми последовательностями"""
    #сколько всего белков создано
    count = 0

    #словарь {id белка : а/к последовательность}
    id_prot = {}

    #контрольная последовательность
    control = "SEGAFGRDFHNRESHJRSTTTTHHREHSEHEHJ"

    def __init__(self, sequence):
        self._sequence = sequence
        Protein.count += 1

        #генерируем уникальный id
        prot_id = random.randint(1, 100)
        while prot_id in Protein.id_prot:
            prot_id = random.randint(1, 100)
        Protein.id_prot[prot_id] = sequence
        self._id = prot_id

    def sequence(self):
        """Возвращает а/к последовательность экземпляра"""
        return self._sequence

    @staticmethod
    def get_proteins():
        """Возвращает словарь всех созданных белков"""
        return Protein.id_prot
    
    @staticmethod
    def reset_proteins():
        """Очищает словарь белков"""
        Protein.id_prot.clear()
        Protein.count = 0

    def homology(self):
        """Возвращает процент гомологичности с контрольной последовательностью"""
        seq = self._sequence.upper()
        ctrl = Protein.control.upper()

        #выравниваем длину под контрольную
        if len(seq) > len(ctrl):
            seq = seq[:len(ctrl)] #обрезаем
        elif len(seq) < len(ctrl):
            seq = seq + "*" * (len(ctrl) - len(seq)) #дополняем

        #считаем гомологичность
        matches = sum(1 for a, b in zip(seq, ctrl) if a == b)
        return matches / len(ctrl) * 100

    def __str__(self):
        return f"Protein(id={self._id}, len={len(self._sequence)})"

    def __repr__(self):
        return f"Protein({self._sequence})"
# %%
### проверки задания 2.5### 
Protein.reset_proteins()   
p1 = Protein("SEGAFGRDFlllESHJRSTTTTHHREHSEHEHJ")
p2 = Protein("SEGAFesfesfgrrhHNRESHJRSTTTTHHREpopopHJrgrdggdg")
p3 = Protein("hhhhh")

print(Protein.count,
Protein.id_prot,
p1.sequence(),
p2.sequence(),
p3.sequence(), sep="\n"
)