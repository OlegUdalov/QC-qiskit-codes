'''this file contains classical subroutines of the Shor algorithm'''

#first is the Euclids algorithm to find greatest common divisor (GCD) of two numbers a and b
def euclids_GCD(a, b):
    '''the function finds the greatest common divisor (GCD) of two numbers a and b'''
    # the input parameters should be positive integers
    found = False
    if (type(a) != int) or (type(b) != int):
        print('Input parameters should be positive integer numebrs!')
        res = -1
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
            
#this function finds the inverse of an integer number b modulo N
def euclids_inverse_mod(a, N):
    '''this function finds the inverse of an integer number b modulo N
    y = a^(-1) mod N
    Parameter:
    __________
    a: int, intger number to be inverted
    N: int, 
    
    
    Returns:
    __________
    int, a^(-1) mod N

    Errors:
    __________
    if function returns -1, then there is an error in the input parameters


    References:
    ___________
    https://en.wikipedia.org/wiki/Modular_multiplicative_inverse
    '''
    
    # the input parameters should be positive integers
    found = False
    print(N)
    print(a)
    if (type(N) != int) or (type(a) != int):
        print('Input parameters should be positive integer numebrs!')
        res = -1
    else:
        if (a - round(a) != 0) or (N - round(N) != 0) or (a <= 0) or (N <= 0):
            print('Input parameters should be positive integer numebrs!')

        r = []
        s = []
        t = []
        q = []
        
        if N > a:
            r.append(N)
            r.append(a)
        else:
            r.append(a)
            r.append(N)
        s.append(1)
        s.append(0)
        t.append(0)
        t.append(1)
        q.append(0)
        
        i = 1
        while r[i] != 0:
            q.append(r[i - 1] // r[i])
            r.append(r[i - 1] - q[i] * r[i])
            s.append(s[i - 1] - q[i] * s[i])
            t.append(t[i - 1] - q[i] * t[i])
            i = i + 1
                  
        if N > a:
            if t[i - 1] < 0:
                t[i - 1] = t[i - 1] + N
            res = t[i-1]
            
        if N < a:
            if s[i - 1] < 0:
                s[i - 1] = s[i - 1] + N
            res =  s[i-1]
    print(res)
    return res

#converting binary number represented by a list of ones and zeroes to a integer
def bin_str_2int(a):
    '''converting binary number represented by a list of ones and zeroes to a integer
    string like [1, 0, 1] is converted to int 5.
    Note that the first element of the list represents the highest digit
    
    Parameter:
    __________
    a: list of zeroes and ones
    
    
    
    Returns:
    __________
    int, binary number converted to int
    '''
    res = 0
    for x in a:
        if x != 0 and x != 1:
            res = -1
            print('Not a binary number')
    if res == 0:
        res = sum(x * 2**i for i, x in enumerate(reversed(a)))
    
    return res

#converting int number to a binary represented by a list of ones and zeroes 
def int_2_bin_str(a,n):
    '''converting int number to a binary represented by a list of ones and zeroes 
    converted to int like 5 to a string like [1, 0, 1]. NUmber of digits in the list is defined by the number n
    Note that the first element of the list represents the highest digit
    
    Parameter:
    __________
    a, int, in number  to convert
    n, int, number of digits we should have

    Errors:
    __________
    1) if the input number is less than 0, then the -1 is returned and text message is given
    
    Warnings:
    1) if the n is less that number of digit n we set, then the warning is provided.
    
    Returns:
    __________
    list, list of zeroes and ones
    '''
    if a < 0:
        b = -1
        print('Error: number should be positive')
    else:
        x = (list(bin(a)))
        del x[0]
        del x[0]
        i = 0
        b = []
        b[:] = x[:]
        for z in x:
            z = int(z)
            b[i] = z
            i = i + 1
        s = len(b)
        if s > n:
            print('Warning: need more digits')
        else:
            if s < n:
                dig = []
                for i in range(n-s):
                    dig.append(0)

                b = dig + b
                
            
       
    return b

