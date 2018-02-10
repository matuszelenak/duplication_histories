#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np
from random import randint
import os
import argparse
import math
import json

i = {'A' : 'T', 'C' : 'G', 'T' : 'A', 'G' : 'C'}

id_counter = 0

def generateID():
	global id_counter
	id_counter += 1
	return id_counter

class Atom():
	def __init__(self, t, id, len):
		self.type = t
		self.id = id
		self.length = len
		self.sequence = np.random.choice(['A', 'C', 'G', 'T'], len, p = [0.25, 0.25, 0.25, 0.25]).tolist()

	def invert(self):
		self.type = - self.type
		self.sequence = [i[b] for b in self.sequence]

	def duplicate(self):
		a = Atom(self.type, generateID(), self.length)
		a.sequence = self.sequence[::]
		return a

	def __str__(self):
		return str(self.type)

	def __repr__(self):
		return str(self.type)

def perform_duplication(atom_seq, centr, l, distance, typ, delet):
	event_params = {}
	event_params['centroid'] = centr
	event_params['length'] = l
	event_params['distance'] = distance
	event_params['type'] = typ
	event_params['direction'] = np.random.choice(['left', 'right'], p = [0.5, 0.5])

	breakpoint = centroid + distance
	l, r = centr - (l - 1)//2 , centr + l // 2 + 1
	duplicated = [a.duplicate() for a in atom_seq[l:r]]
	if typ == 'dup_inv':
		for x in duplicated:
			x.invert()
	if delet and len(duplicated) > 2:
		c2 = randint(0, len(duplicated) - 1)
		len2 = len(duplicated) + 1
		while not ((len2 - 1)//2 <= c2 and len2 // 2 <= len(duplicated) - c2 - 1):
			len2 = np.random.geometric(1/10)
		l2, r2 = c2 - (len2 - 1)//2 , c2 + len2 // 2 + 1
		duplicated = duplicated[:l2] + duplicated[r2:]

		deletion = {}
		deletion['centroid'] = c2
		deletion['length'] = len2
		event_params['deletion'] = deletion

	if event_params['direction'] == 'left':
		duplicated = duplicated[::-1]

	return atom_seq[:breakpoint] + duplicated + atom_seq[breakpoint:], event_params


parser = argparse.ArgumentParser()

args = parser.parse_args()

#mean_length = 14307
#mean_distance = 306718
p_dup_inv = 0.39
p_dup = 1 - p_dup_inv
p_del = 0.1

for test_sample in range(1000):
	j = {}
	events = []
	id_counter = 0
	num_of_atoms = 	randint(5, 30)
	atom_lengths = [randint(20, 1000) for _ in range(num_of_atoms)]
	#print(atom_lengths)
	ancestral_sequence = [Atom(t, generateID(), atom_lengths[t]) for t in np.random.permutation(num_of_atoms)]
	num_of_events = randint(1, 30)
	events_types = np.random.choice(['dup', 'dup_inv'], num_of_events, p = [p_dup, p_dup_inv])
	events_dels = np.random.choice(['del', 'no_del'], num_of_events, p = [p_del, 1 - p_del])

	atom_sequence = ancestral_sequence
	for e, d in zip(events_types, events_dels):
		centroid = randint(0, len(atom_sequence) - 1)
		length = len(atom_sequence) + 1
		while not ((length - 1)//2 <= centroid and length // 2 <= len(atom_sequence) - centroid - 1):
			length = np.random.geometric(1/10)

		l, r = centroid - (length - 1)//2 , centroid + length // 2 + 1
		distance = len(atom_sequence) + 1
		while True:
			distance = np.random.geometric(1/20) - 1
			if (r + distance <= len(atom_sequence)):
				distance += length // 2 + 1
				break
			if (l - distance >= 0):
				distance = - distance - (length - 1)//2
				break
		atom_sequence, event_params = perform_duplication(atom_sequence, centroid, length, distance, e, d == 'del')
		events.append(event_params)
	j['events'] = events
	j['ancestral'] = [{"type" : int(x.type), "sequence" : "".join([y for y in x.sequence])} for x in ancestral_sequence]
	j['contemporary'] = [{"type" : int(x.type), "sequence" : "".join([y for y in x.sequence])} for x in atom_sequence]
	f = open(os.path.join("./simulated_data", 'history_{0}.json'.format(test_sample)), 'w')
	f.write(json.dumps(j, indent=4, sort_keys=True, ensure_ascii=False))
