open('Contacts-6hrs.txt', 'w').close()
contacts = open("Contacts-6hrs.txt","a")
idList = [line.rstrip('\r\n') for line in open("contact_ids.txt")]
print(idList)

for name_id in idList:
	genes = open("genes6.txt","r")
	for line in genes:
		data_list = line.split()
		if name_id == data_list[1]:
			contacts.write(line)
	