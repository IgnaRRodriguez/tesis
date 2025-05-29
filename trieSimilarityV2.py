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


def addProtein(node:Node,protein:str,family:str):
	currentNode = node
	level = 0
	for letter in protein:
		position= ord(letter) - 65
		if currentNode.nexts[position] is None :
			currentNode.nexts[position]= Node()
		currentNode = currentNode.nexts[position]
		if family not in currentNode.families: 
			currentNode.families.add(family)
		level += 1

def getDistanceFamily(trie:Node,family1:str,family2:str,totalNodes:dict):
	nextTries = deque()
	similarty = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families) and (family2 in nexts.families)):
			nextTries.append(nexts)
	level = 0
	while nextTries:
		currentTries = nextTries 
		nextTries = deque()
		level += 1
		levelSimilarity = 0
		while currentTries:
			currentNode = currentTries.popleft()
			nodeScore = 1
			for family in currentNode.families:
				nodeScore *= (totalNodes[level][family]/totalNodes[level]['total'])
			levelSimilarity += nodeScore
			for nexts in filter(None,currentNode.nexts):
				if ((family1 in nexts.families) and (family2 in nexts.families)):
					nextTries.append(nexts)
		similarty += levelSimilarity
	return similarty
	


def pruneTrie(trie:Node,totalNodes:dict):
	level = 0
	for i in range(len(trie.nexts)):
		if (trie.nexts[i] is not None and len(trie.nexts[i].families) ==1 ):
			trie.nexts[i] = None

	nextLevel = deque(trie.nexts)
	while nextLevel:
		level += 1
		totalNodes[level] = {}
		totalNodes[level]['total'] = 0
		currentLevel = nextLevel
		nextLevel= deque()
		while currentLevel:
			currentNode = currentLevel.popleft()
			if currentNode is not None:
				totalNodes[level]['total'] += 1
				for i in range(len(currentNode.nexts)):
					if (currentNode.nexts[i] is not None and len(currentNode.nexts[i].families) == 1 ):
						currentNode.nexts[i] = None
				for family in currentNode.families:
					if family in totalNodes[level]:
						totalNodes[level][family] += 1
					else:
						totalNodes[level][family] = 1
				nextLevel.extend(currentNode.nexts)


print("start")
ignoredFamilies = {}
trie= getTrie()
dirList = list(filter(lambda x: (x.endswith("NE.csv")),os.listdir("../trie/output")))
families = []
totalNodesByLevel = {}
proteinsInFamily = {} 
with open('../familiesInfo.csv',mode='r') as csv_file:
	csv_reader = csv.DictReader(csv_file)
	for row in csv_reader:
		proteinsInFamily[row['family']]=int(row['quantity'])
startTime=datetime.datetime.now()
print("get trie")
for family in dirList:
	with open('../trie/output/'+family,mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		name= family.split("_")[0]+"_"+family.split("_")[1]
		if name not in ignoredFamilies:
			families.append(name)

			for row in csv_reader:
				#Pattern,PatternLength,Instances,AffectedProteins
				if ((int(row["PatternLength"])<=25) and ((int(row["AffectedProteins"])/proteinsInFamily[name])>=0.05)):
					addProtein(trie,row["Pattern"],name)

pruneTrie(trie,totalNodesByLevel)
print(totalNodesByLevel)

distances= [([0]* len(families)) for i in range(len(families))]
for i in range(len(families)):
	family1 = families[i]
	print(family1)
	for j in range(i,len(families)):
		family2 = families[j]
		print(family2)
		val = getDistanceFamily(trie,family1,family2,totalNodesByLevel)
		distances[i][j] = val
		distances[j][i] = distances[i][j]
	print("-------------------------------------------------")


endTime=datetime.datetime.now()
print('time: ',abs(endTime-startTime).total_seconds() / 60.0)
with open('resultV3.2-ProteinasAfectadas-05Short.csv','w',newline='') as csv_file:
	csvWriter = csv.writer(csv_file)
	csvWriter.writerow(families)
	csvWriter.writerows(distances)