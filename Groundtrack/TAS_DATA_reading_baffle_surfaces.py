#IMPORT MODULES
from matplotlib import pyplot as plt

#READ FILE AND CREATE LIST WITH DATA
r = open('baffle_surfaces.csv','r',encoding='utf-8',  errors='replace')
data = []
for line in r.readlines():  #saves all data each line becomes a list
    line = line.strip()
    line = line.strip('\n')
    lst  = line.split(',')
    data.append(lst)
r.close()

#TRANSFORM DATA LIST INTO FLOATS
for i in range(6, len(data)):
   for j in range(0, len(data[i])):
       data[i][j] = float(data[i][j])

'''   THIS IS JUST A CHECK POINT 
#CHECK - data from 6 to 807
print(data[6])  #first data line
print(data[len(data)-1])  #last data line
print(type(data[6][1]))  #now is a float!!!!!
print(type(data[len(data)-1][1]))
'''

#GET NAMES OF ELEMENTS OF TELESCOPE FROM LINE 2 OF LIST CREATED
names = []
for i in range(1,len(data[2]),2):
    names.append(data[2][i])
print(names)   #this is a list with the names of the 8 data locations,
# for each data location there is a time and a temperature recorded

#FORGET ALL THE LINES THAT DO NOT CONTAIN DATA (BEFORE LINE 6 IS USELESS NOW)
data = data[6:]

#CREATE A LIST FOR EACH ELEMENT WITH STRUCTRE: [[NAME], [X1,Y1],[X2,Y2],...]
#First create list with name already input
lst0 = [names[0]]
lst1 = [names[1]]
lst2 = [names[2]]
lst3 = [names[3]]
lst4 = [names[4]]
lst5 = [names[5]]
lst6 = [names[6]]
lst7 = [names[7]]

#Assign the relevant columns (x,y) to the elements  [check top lines of excel file for clarification]
for line in data:
     lst0.append([line[0],line[1]])
     lst1.append([line[2],line[3]])
     lst2.append([line[4],line[5]])
     lst3.append([line[6],line[7]])
     lst4.append([line[8],line[9]])
     lst5.append([line[10],line[11]])
     lst6.append([line[12],line[13]])
     lst7.append([line[14],line[15]])

#These are the Time-Temp measurements for each element
print(lst0)
print(lst1)
print(lst2)
print(lst3)
print(lst4)
print(lst5)
print(lst6)
print(lst7)


#DEFINE FUNCTION TO PRINT TEMP(TIME) FOR EACH ELEMENT
def plot(list):
    #coordinate lsts
    x = []
    y = []
    for i in range(1,len(list)):  #given that structure of each list is [[NAME], [X1,Y1],[X2,Y2],...]
        x.append(list[int(i)][0])
        y.append(list[int(i)][1])

    #plot
    plt.figure(figsize=(14, 10))
    plt.plot(x,y, 'x')
    plt.title('Element  ' + list[0]) #the title is taken from the item 0 of the list which is name of element
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [C]')
    plt.grid()

    plt.show()

#PLOT Data taking: plot(list, list[0](which is the title))
plot(lst0)
plot(lst1)
plot(lst2)
plot(lst3)
plot(lst4)
plot(lst5)
plot(lst6)
plot(lst7)
#I just plotted the first two lists as an example

