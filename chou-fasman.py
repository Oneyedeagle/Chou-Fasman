import numpy as np

#function for assigning propensity
def assign_propensity(residue, protein_type):
	propensity=0
	if (residue == 'A'):
		propensity = np.array([1.45, 0.97])
	elif (residue == 'C'):
		propensity = np.array([0.77, 1.30])
	elif (residue == 'D'):
		propensity = np.array([0.98, 0.80])
	elif (residue == 'E'):
		propensity = np.array([1.53, 0.26])
	elif (residue == 'F'):
		propensity = np.array([1.12, 1.28])
	elif (residue == 'G'):
		propensity = np.array([0.53, 0.81])
	elif (residue == 'H'):
		propensity = np.array([1.24, 0.71])
	elif (residue == 'I'):
		propensity = np.array([1.00, 1.60])
	elif (residue == 'K'):
		propensity = np.array([1.07, 0.74])
	elif (residue == 'L'):
		propensity = np.array([1.34, 1.22])
	elif (residue == 'M'):
		propensity = np.array([1.20, 1.67])
	elif (residue == 'N'):
		propensity = np.array([0.73, 0.65])
	elif (residue == 'P'):
		propensity = np.array([0.59, 0.62])
	elif (residue == 'Q'):
		propensity = np.array([1.17, 1.23])
	elif (residue == 'R'):
		propensity = np.array([0.79, 0.90])
	elif (residue == 'S'):
		propensity = np.array([0.79, 0.72])
	elif (residue == 'T'):
		propensity = np.array([0.82, 1.20])
	elif (residue == 'V'):
		propensity = np.array([1.14, 1.65])
	elif (residue == 'W'):
		propensity = np.array([1.14, 1.19])
	elif (residue == 'Y'):
		propensity = np.array([0.61, 1.29])

	if(protein_type == 'helix'):                             #0th index of array has propensity of helices
		return (propensity[0])

	elif(protein_type == 'beta'):                            #1st index of array has propensity of strand
		return (propensity[1])


#------------------------------------------------

#e.g when i = 4 and protein is SGFRKMAFPSG
#  0 1 2 3 4 5 6 7 8 9 10
# -----------------------
# |S|G|F|R|K|M|A|F|P|S|G|
# ----------------------- 
#          ^   ^ ^
#          i   k l

#------------------------------------------------


#extend for alpha helix
def backward(x, i):
	while(i>=4):
		if(assign_propensity(x[i-1], 'helix') + assign_propensity(x[i-2], 'helix') + assign_propensity(x[i-3], 'helix') + assign_propensity(x[i-4], 'helix')>=4):
			i = i-1					#if extend == true; assign l the value of k = (i+2) -1 = i
		else:
			return i                #returns when Σpropensity is not >=4
	return i


def forward(x, i):
	while(i<=len(x)-5):
		if(assign_propensity(x[i+1], 'helix') + assign_propensity(x[i+2], 'helix') + assign_propensity(x[i+3], 'helix') + assign_propensity(x[i+4], 'helix')>=4):
			i = i+1					#if extend == true; assign k the value of k = (i+2) +1 = i+3
		else:
			return i  				#returns when Σpropensity is not >=4
	return i

#extend for beta strand
def backward2(x, i):
	while(i>=4):
		if(assign_propensity(x[i-1], 'beta') + assign_propensity(x[i-2], 'beta') + assign_propensity(x[i-3], 'beta') + assign_propensity(x[i-4], 'beta')>=4):
			i = i-1
		else:
			return i
	return i


def forward2(x, i):
	while(i<=len(x)-5):
		if(assign_propensity(x[i+1], 'beta') + assign_propensity(x[i+2], 'beta') + assign_propensity(x[i+3], 'beta') + assign_propensity(x[i+4], 'beta')>=4):
			i = i+1					#if extend = true assign k the value of k = (i+2) +1 = i+3
		else:
			return i
	return i

#to find the 6 residues who have Σpropensity >=4
def find_helix(x,N):

	i=0						                                #starting index of 6 residue
	l=0						                                #index for backward extension
	k=0						                                #index for forward extension
	alphaans=[]                                             #to store alpha list
	for j in range(0, len(x)):
		alphaans.append(" ")
	score = 0
	while (i<N-5):
		count = 0                                           #counter
		for j in range(i,i+6):								#summation of 6 window (residues)
			if(assign_propensity(x[j], 'helix')>=1):
				count = count +1

		if(count>=4.0):										#if true: 6 residues found who have score >=4
			# print(x[i:i+6])
			flag1=0
			flag2=0
			l=i+3											#starting index of left extension
			k=i+2                                           #ending index of right extension
			
			if(i>=4):
				l = backward(x,l)							#finding backward extension's index

			if(i<N-5):
				k = forward(x,k)							#finding forward extension's index

			if(l!=i+3):										#if true implies extended backward(left)
				# print(l-3,end=' ')						#printing the starting index of left extended helix
				flag1=1
				l = l-3
			else:
				# print(i, end=' ')							#again printing the starting index of left extended helix if it doesnt get extended
				pass

			if(k!=i+2):
				flag2=1
				# print(k+3)								#printing the ending index of right extended helix
				k = k+3										#k is assigned the ending index
			else:
				# print(i+5)								#again printing the starting index of right extended helix if it doesnt get extended
				pass

			if flag1==0:
				l = i
			if flag2 ==0:
				k=i+5
			# print(x[l:k+1])								#printing the whole helix

			for j in range(l,k+1):							
				alphaans[j]='H'



		i = i+1
	return alphaans

#to find the 5 residues who have Σpropensity >=4
def find_strand(x,N):

	i=0	                                                #starting index of 5 residue
	l=0						                            #index for backward extension
	k=0						                            #index for forward extension
	strandans=[]                                         #to store beta list
	for j in range(0, len(x)):
		strandans.append(" ")
	score = 0
	while (i<N-4):
		store = 0
		count = 0
		for j in range(i,i+4):							#Σpropensity of 5 window (residues)
			if(assign_propensity(x[j], 'beta')>=1):
				count = count +1
	
		if(count>=3):									#if true: 6 residues found who have score >=4
			
			flag1=0
			flag2=0
			l=i+3										#ending index __ of left extension
			k=i+1
			if(i>=4):

				l = backward2(x,l)						#finding backward extension's index

			if(i<=N-5):

				k = forward2(x,k)						#finding forward extension's index

			if(l!=i+3):									#if true implies extended backward(left)
				# print(l-3,end=' ')					#printing the starting index of left extended helix
				flag1=1									#flag is ON if exxtended == true
				l = l-3
			else:
				# print(i, end=' ')						#again printing the starting index of left extended helix if it doesnt get extended
				pass

			if(k!=i+1):
				flag2=1									#flag2 is ON if right extended == true
				# print(k+3)							#printing the ending index of right extended helix
				k = k+3									#k is assigned the ending index
			else:
				# print(i+5)							#again printing the starting index of right extended helix if it doesnt get extended
				pass

			if flag1==0:
				l = i
			if flag2 ==0:
				k=i+4
			# print(x[l:k+1])							#printing the whole strand
			for j in range(l,k+1):
				strandans[j]='S'



		i = i+1
	return strandans
	
#to convert list -> string
def convert(ans):
	y = ""
	for i in ans:
		y += i
	return y
		


def main():

	x = str(input("Enter the Protein sequence: "))
	N = len(x)
	

	ans1 = find_helix(x,N)																	#Helix list i.e ['H', ' ', 'H', 'H' , 'H'...]
	ans2 = find_strand(x,N)																	#Beta list i.e ['S', ' ', 'S', 'S' , 'S'...]

	sum_alpha = 0
	sum_beta = 0


	ans3 = []
	
	i=0
	while (i <len(x)):
		if(ans1[i] != ' '):																	#if helix list's ith pos. is not empty					
			if(ans2[i] != ' '):																#if beta list's ith pos. is not empty			
				
				start_val = i 																#stores the starting index of conflicting sequences
				end_val = i  																#stores the ending index of conflicting sequences 
				
				while((ans1[i] != ' ' ) and (ans2[i] != ' ')): 								#loop condn where we have conflicting sequences (helix or strand?)
					sum_alpha = sum_alpha + assign_propensity(x[i], 'helix')
					sum_beta = sum_beta + assign_propensity(x[i], 'beta')
					end_val = end_val+1
					i = i+1

				if(sum_alpha>=sum_beta):													#if Σpropensity (helix) > Σpropensity (beta)
					for j in range(start_val, end_val):
						ans3.append('H')
						sum_alpha = 0  														#resetting the value of sum_alpha
						sum_beta = 0  														#resetting the value of sum_beta


				elif(sum_alpha<sum_beta):                                                   #if Σpropensity (helix) < Σpropensity (beta)
					for j in range(start_val, end_val):                                    
						ans3.append('S')
						sum_alpha = 0														#resetting the value of sum_alpha
						sum_beta = 0  														#resetting the value of sum_beta


			elif(ans2[i] == ' '): 															#if list1 has alpha helix and list2 is empty
				ans3.append(ans1[i])  														#the reslatant list will have alpha helix at that position
				i = i+1


		elif(ans1[i] == ' '):  																#if list1 is empty and list2 has beta strand at ith pos
			if(ans2[i] != ' '):      														#the reslatant list will have beta beta at that position
				ans3.append(ans2[i])
				i = i+1

			elif(ans2[i] == ' '):                                                           #if both lists are empty at ith pos
				ans3.append('T')                                                            #resultant list will have a Turn there
				i = i+1
	

	ans_1 = convert(ans1)							                                        #Helix string i.e H HHHHH...
	ans_2 = convert(ans2)							                                        #Beta string i.e S SSSSS...		
	ans_3 = convert(ans3)							                                        #Resultant string


#printing the 4 strings
	# print(x)											                                    #the original protein
	print("                           ",ans_1)                                                                            #Helix strand
	print("                           ",ans_2)                                                                            #Beta strand
	print("                           ",ans_3)                                                                            #resultant strand
	

if __name__ == '__main__':
	main()
