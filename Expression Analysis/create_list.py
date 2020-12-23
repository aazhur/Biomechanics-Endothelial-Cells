open('allgenes.txt', 'w').close()
contacts = open("allgenes.txt","a")

genes = open("genes6.txt","r")
for line in genes:
	data_list = line.split()
	contacts.write(data_list[2])
	contacts.write('\r\n')
			
contacts.close()		
idList = [line.rstrip('\r\n') for line in open("allgenes.txt","r")]	
contacts = open("allgenes.txt","a")		

genes = open("genes2.txt","r")
for line in genes:
	data_list = line.split()
	score = 0
	for name_id in idList:
		if name_id == data_list[2]:
			score += 1	
	if score == 0:		
		contacts.write(data_list[2])
		contacts.write('\r\n')
		
contacts.close()			
idList = [line.rstrip('\r\n') for line in open("allgenes.txt","r")]	
contacts = open("allgenes.txt","a")		

genes = open("genes0.txt","r")
for line in genes:
	data_list = line.split()
	score = 0
	for name_id in idList:
		if name_id == data_list[2]:
			score += 1	
	if score == 0:		
		contacts.write(data_list[2])
		contacts.write('\r\n')	

contacts.close()			