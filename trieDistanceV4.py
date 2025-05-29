import csv
import os
import datetime
import math
import json
from collections import deque

class Node(object):
	"""docstring for Node"""
	def __init__(self,):
		self.families = {}
		self.nexts= [None] * 26
		
def getTrie():
	return Node()


def addProtein(node:Node,protein:str,family:str,affectedProtein:int):
	currentNode = node
	level = 0
	for letter in protein:
		position= ord(letter) - 65
		if currentNode.nexts[position] is None :
			currentNode.nexts[position]= Node()
		currentNode = currentNode.nexts[position]
		if family not in currentNode.families.keys(): 
			currentNode.families[family]=0
		level += 1
	currentNode.families[family]=affectedProtein

def getDistanceFamily(trie:Node,family1:str,family2:str,familiesInfo:dict):
	nextTries = deque()
	distance = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families.keys()) or (family2 in nexts.families.keys())):
			nextTries.append(nexts)
	level = 0
	while nextTries:
		currentTries = nextTries 
		nextTries = deque()
		level += 1
		while currentTries:
			currentNode = currentTries.popleft()
			multiplier = 1
			#multiplier = 1 if len(currentNode.families.keys()) == 1 else -math.log(1-(len(currentNode.families.keys())/(len(familiesInfo.keys())+1)),2) 
			if((family1 in currentNode.families.keys()) and (family2 in currentNode.families.keys())):
				family1Frequency = (familiesInfo[family1]['length']+(familiesInfo[family1]['quantity']*(1-level)))
				family2Frequency = (familiesInfo[family2]['length']+(familiesInfo[family2]['quantity']*(1-level)))
				distance += (abs((currentNode.families[family1]/family1Frequency)-(currentNode.families[family2]/family2Frequency))*multiplier)
			elif (family1 in currentNode.families.keys()):
				family1Frequency = (familiesInfo[family1]['length']+(familiesInfo[family1]['quantity']*(1-level)))
				distance += ((currentNode.families[family1]/family1Frequency)*multiplier)
			else:
				family2Frequency = (familiesInfo[family2]['length']+(familiesInfo[family2]['quantity']*(1-level)))
				distance += ((currentNode.families[family2]/family2Frequency)*multiplier)
			for nexts in filter(None,currentNode.nexts):
				if ((family1 in nexts.families.keys()) or (family2 in nexts.families.keys())):
					nextTries.append(nexts)
		
	return distance

def pruneTrie(trie:Node,totalNodes:dict):
	level = 0
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

startTime=datetime.datetime.now()

with open('familiesInfo.json', 'r') as file:
    familiesInfo = json.load(file)


print("get trie")
for family in dirList:
	with open('../trie/output/'+family,mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		name= family.split("_")[0]+"_"+family.split("_")[1]
		if name not in ignoredFamilies:
			families.append(name)
			for row in csv_reader:
				#Pattern,PatternLength,Instances,AffectedProteins <=10
				if ((int(row["PatternLength"])<6)):
					addProtein(trie,row["Pattern"],name,int(row["Instances"]))
#totalNodesByLevel={}
#pruneTrie(trie,totalNodesByLevel)
#print(totalNodesByLevel)


distances= [([0]* len(families)) for i in range(len(families))]
for i in range(len(families)):
	family1 = families[i]
	print(family1)
	for j in range(i,len(families)):
		family2 = families[j]
		print(family2)
		val = getDistanceFamily(trie,family1,family2,familiesInfo)
		distances[i][j] = val
		distances[j][i] = distances[i][j]
	print("-------------------------------------------------")


endTime=datetime.datetime.now()
print('time: ',abs(endTime-startTime).total_seconds() / 60.0)
with open('resultV4.csv','w',newline='') as csv_file:
	csvWriter = csv.writer(csv_file)
	csvWriter.writerow(families)
	csvWriter.writerows(distances)