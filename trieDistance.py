import csv
import os
import datetime
import math
from collections import deque

class Node(object):
	"""docstring for Node"""
	def __init__(self,):
		self.families = set()
		self.nexts= [None] * 26
		
def getTrie():
	return Node()


def addProtein(node:Node,protein:str,family:str,totalNodes:dict):
	currentNode = node
	for letter in protein:
		position= ord(letter) - 65
		if currentNode.nexts[position] is None :
			totalNodes['total'] += 1
			currentNode.nexts[position]= Node()
		currentNode = currentNode.nexts[position]
		if family not in currentNode.families: 
			currentNode.families.add(family)
			totalNodes[family] += 1

def getDistanceFamily(trie:Node,family1:str,family2:str,totalNodes:dict):
	nextTries = deque()
	diff = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families) or (family2 in nexts.families)):
			nextTries.append(nexts)

	while nextTries:
		nextTrie = nextTries.popleft()
		if ((family1 in nextTrie.families) and (family2 in nextTrie.families)):
			diff += (1- (len({family1,family2})/len(nextTrie.families)))
		else:
			diff += 1
		for nexts in filter(None,nextTrie.nexts):
			if ((family1 in nexts.families) or (family2 in nexts.families)):
				nextTries.append(nexts)
	div = totalNodes[family1] + totalNodes[family2]
	return diff/div	

def getDistanceFamilyNew(trie:Node,family1:str,family2:str,totalNodes:dict):
	nextTries = deque()
	diff = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families) or (family2 in nexts.families)):
			nextTries.append(nexts)

	while nextTries:
		nextTrie = nextTries.popleft()
		if ((family1 in nextTrie.families) and (family2 in nextTrie.families)):
			diff += abs(distance1(nextTrie,totalNodes,family1)-distance1(nextTrie,totalNodes,family2))
			#diff += abs(distance2(nextTrie,totalNodes,family1)-distance2(nextTrie,totalNodes,family2))
		else:
			diff += 1
		for nexts in filter(None,nextTrie.nexts):
			if ((family1 in nexts.families) or (family2 in nexts.families)):
				nextTries.append(nexts)
	div = totalNodes[family1] + totalNodes[family2]
	return diff/div

def getDistanceFamilyNew2(trie:Node,family1:str,family2:str,totalNodes:dict):
	nextTries = deque()
	diff = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families) or (family2 in nexts.families)):
			nextTries.append(nexts)

	while nextTries:
		nextTrie = nextTries.popleft()
		if ((family1 in nextTrie.families) or (family2 in nextTrie.families)):
			diff += abs(distance3(nextTrie,totalNodes,family1)-distance3(nextTrie,totalNodes,family2))
			#diff += abs(distance4(nextTrie,totalNodes,family1)-distance4(nextTrie,totalNodes,family2))
		for nexts in filter(None,nextTrie.nexts):
			if ((family1 in nexts.families) or (family2 in nexts.families)):
				nextTries.append(nexts)
	div = totalNodes[family1] + totalNodes[family2]
	print(diff)
	return diff/div	

def distance1(node:Node,totalNodes:dict,family:str):
	prob = totalNodes[family]/totalNodes['total']
	cantFams = len(node.families)
	base = cantFams
	if(cantFams == 1):
		return 0
	else:
		result = ((-1)*cantFams*(prob*math.log(prob,base)))
		return result

def distance2(node:Node,totalNodes:dict,family:str):
	prob = totalNodes[family]/totalNodes['total']
	cantFams = len(node.families)
	base = totalNodes['total']
	if(cantFams == 1):
		return 0
	else:	
		return ((-1)*cantFams*(prob*math.log(prob,base)))

def distance3(node:Node,totalNodes:dict,family:str):
	prob = totalNodes[family]/totalNodes['total']
	cantFams = len(node.families)
	base = cantFams
	if(cantFams == 1):
		return int( not (family in node.families))
	elif not (family in node.families):
		return 0
	else:	
		result = ((-1)*cantFams*(prob*math.log(prob,base)))
		return result

def distance4(node:Node,totalNodes:dict,family:str):
	prob = totalNodes[family]/totalNodes['total']
	base = totalNodes['total']
	cantFams = len(node.families)
	if(cantFams == 1):
		return int(not(family in node.families))
	else:	
		return ((-1)*cantFams*(prob*math.log(prob,base)))

print("start")
trie= getTrie()
dirList = list(filter(lambda x: x.endswith("NE.csv"),os.listdir("../trie/output")))
families = []
totalNodes = {'total':0}
startTime=datetime.datetime.now()
print("get trie")
for family in dirList:
	with open('../trie/output/'+family,mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		name= family.split("_")[0]+"_"+family.split("_")[1]
		families.append(name)
		totalNodes[name] = 0
		for row in csv_reader:
			#Pattern,PatternLength,Instances,AffectedProteins
			if int(row["PatternLength"])<15:
				addProtein(trie,row["Pattern"],name,totalNodes)

distances= [([0]* len(families)) for i in range(len(families))]
print(totalNodes)
for i in range(len(families)):
	family1 = families[i]
	print(family1)
	for j in range(i,len(families)):
		family2 = families[j]
		print(family2)
		val = getDistanceFamilyNew2(trie,family1,family2,totalNodes)
		print(val)
		distances[i][j] = val
		distances[j][i] = distances[i][j]
	print("-------------------------------------------------")


endTime=datetime.datetime.now()
print('time: ',abs(endTime-startTime).total_seconds() / 60.0)
with open('result.csv','w',newline='') as csv_file:
	csvWriter = csv.writer(csv_file)
	csvWriter.writerow(families)
	csvWriter.writerows(distances)