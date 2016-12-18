import gmpy2
from time import time


def new_prime(n_bits=1024):
    """
    is_prime: Miller-Rabin's test (default=25 times)
    """
    upper_bound = gmpy2.mpz(2**n_bits)

    state = gmpy2.random_state(int_time())
    rand_num = gmpy2.mpz_random(state,upper_bound)
    rand_num = gmpy2.bit_set(rand_num,n_bits-1)
    while not gmpy2.is_prime(rand_num):
        state = gmpy2.random_state(int_time())
        rand_num = gmpy2.mpz_random(state,upper_bound)
        rand_num = gmpy2.bit_set(rand_num,n_bits-1)
    return rand_num



def int_time():
    return int(round(time() * 1000))

class Paillier(object):
    @staticmethod
    def generate_keys(n_bits=1024):
        """
        Generates keys necessary to Pailler's Cryptossystem 
        
        Parameters:
            Kwargs:
                n_bits: Size of the key, in bits (default=1024)
        
        Returns:
            Tuple: Tuple containing two tuples:
                t0: Tuple containing the public key (n,g)
                t1: Tuple containing the private key (lambda,mu)
                   
        """
        
        one = gmpy2.mpz(1)
        p = new_prime(n_bits=n_bits)
        q = new_prime(n_bits=n_bits)
        
        # public key
        n = gmpy2.mul(p,q)
        g = gmpy2.add(n,one)
        
        #private key
        lbd = gmpy2.mul(gmpy2.sub(p,one),gmpy2.sub(q,one))
        mu = gmpy2.powmod(lbd,-1,n)
        
        
        return ((n,g),(lbd,mu))
    
    @staticmethod
    def encrypt(m,n,g):
        """
        Encrypts a message, given a public key.
        
        Parameters:
            m: number to be encrypted
            n,g: public keys used for the decryption
        Returns:
            c: encrypted representation of m
        """
        n2 = pow(n,2)
        one = gmpy2.mpz(1)
        state = gmpy2.random_state(int_time())
        r = gmpy2.mpz_random(state,n)
        while gmpy2.gcd(r,n) != one:
            state = gmpy2.random_state(int_time())
            r = gmpy2.mpz_random(state,n)
        
        x = gmpy2.powmod(r,n,n2)
        c = gmpy2.f_mod(gmpy2.mul(gmpy2.powmod(g,m,n2),x),n2)
        
        return c
        
    @staticmethod
    def decrypt(c,lbd,mu,n):
        """
        Decrypts encrypted number c with respect to private keys
                 lambda and mu (and public key n)
        
        Parameters:
            m: number to be encrypted
            lbd,mu: private keys used for the decryption
            n: public key necessary for decryption
        Returns:
            c: encrypted representation of m
        """
        n2 = pow(n,2)
        one = gmpy2.mpz(1)
        x = gmpy2.sub(gmpy2.powmod(c,lbd,n2),one)
        m = gmpy2.f_mod(gmpy2.mul(gmpy2.f_div(x,n),mu),n)
        
        if m >= gmpy2.f_div(n,2):
            m = m - n
        return m
    
    @staticmethod
    def add(c0,c1):
        """
        Adds two encrypted messages together
        
        Parameters:
            c0,c1: encripted messages to be added together
        Returns:
            c: encrypted sum of both messages
            
        """
        
        return gmpy2.mul(c0,c1)
    @staticmethod
    def encrypted_sum(V):
        """
        Parameters:
            V: vector of encrypted messages to be added up
        Returns:
            Encrypted sum of the vector
        
        """
        
        return reduce(Pailler.add,V)
        
        
    @staticmethod
    def sub(c0,c1,n):
        inv_c1 = gmpy2.invert(c1,n**2)
        return gmpy2.mul(c0,inv_c1)
    @staticmethod
    def mul_plain(c0,m0):
        return pow(c0,m0)
        
