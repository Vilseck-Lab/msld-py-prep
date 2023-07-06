with open("CGenFF_atomtypes.txt",'r') as f:
    current_ats = f.readlines()

current_ats = [at.strip() for at in current_ats]
print(current_ats)
