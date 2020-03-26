#IMPORT MODULES
from matplotlib import pyplot as plt

def choose_dataset(filename):
    #READ FILE AND CREATE LIST WITH DATA
    r = open(filename,'r',encoding='utf-8',  errors='replace')
    data = []
    for line in r.readlines():  #saves all data each line becomes a list
        line = line.strip().strip('\n')
        lst  = line.split(',')
        data.append(lst)
    r.close()

    #TRANSFORM DATA LIST INTO FLOATS
    for i in range(6, len(data)):
        for j in range(0, len(data[i])):
            data[i][j] = float(data[i][j])

    #GET NAMES OF ELEMENTS OF TELESCOPE FROM LINE 2 OF LIST CREATED
    names = []
    for i in range(1,len(data[2]),2):
        names.append(data[2][i])
    #print(names)   #this is a list with the names of the 8 data locations,
    # for each data location there is a time and a temperature recorded

    #FORGET ALL THE LINES THAT DO NOT CONTAIN DATA (BEFORE LINE 6 IS USELESS NOW)
    data = data[6:]
    num_datasets = int(len(data[0])/2) # Calculates the number of datasets within file (number of elements for which there is temp data)

    #CREATE A LIST FOR EACH ELEMENT WITH STRUCTRE: [[NAME], [X1,Y1],[X2,Y2],...]
    #First create list with name already input

    data_tab = [] #1st element correspond to the old lst0, etc...
    for i in range(num_datasets):
        data_tab.append([names[i]])

    #Assign the relevant columns (x,y) to the elements  [check top lines of excel file for clarification]
    for line in data:
        for i in range(len(data_tab)):
            data_tab[i].append([line[i*2],line[i*2 + 1]])

    return data_tab