'''this file contains clasical subroutines of the Shor algorithm'''

#first is the Euclids algorithm to find greatest common divisor (GCD) of two numbers a and b
def euclids_GCD(a,b):
    '''the function finds the greatest common divisor (GCD) of two numbers a and b'''
    # the input parameters should be positive integers
    found = False
    if (type(a)!=int) or (type(b)!=int):
        print('Input parameters should be positive integer numebrs!')
        res=-1
    else:
        if (a - round(a) != 0) or (b - round(b) != 0) or (a <= 0) or (b <= 0):
            print('Input parameters should be positive integer numebrs!')
            res = -1
            found = True
        while found == False:
            if a > b:
                a = a - b
                if a == 0:
                    found = True
                    res = b
                else:
                   if a == 1:
                       found = True
                       res = 1
            else:
                b = b-a
                if b == 0:
                    found = True
                    res = a
                else:
                   if b == 1:
                       found = True
                       res = 1
    return res
            

