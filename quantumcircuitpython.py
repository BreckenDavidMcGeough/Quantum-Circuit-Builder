import numpy as np 


class QCircuit:
    def __init__(self,basis_state_one,basis_state_two,bool_func,dim=2):
        self.basis_state_one = basis_state_one
        self.basis_state_two = basis_state_two
        self.bool_func = bool_func
        self.hadamard_constant = (1/(2**(1/2)))
        self.hadamard_gate = np.array([[1,1],[1,-1]])
        self.I = np.array([[1,0],[0,1]])
        self.NOT_gate = np.array([[0,1],[1,0]])
        self.CNOT_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        self.SWAP_gate = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
        self.AT_gate = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
        
    def deutsch_algorithm(self):
        basis_tensor = np.kron(self.basis_state_one,self.basis_state_two)
        I_tensor_X = np.kron(self.I,self.NOT_gate)
        a1 = np.dot(I_tensor_X,basis_tensor)
        H_tensor_H = np.kron(self.hadamard_gate,self.hadamard_gate)
        self.hadamard_constant = 1/2
        a2 = np.dot(H_tensor_H,a1)
        a3 = np.dot(self.bool_func,a2)
        H_tensor_I = np.kron(self.hadamard_gate,self.I)
        self.hadamard_constant = 1/(2*(2**(1/2)))
        circuit_output = self.hadamard_constant*np.dot(H_tensor_I,a3)
        return circuit_output

    def qbit_entanglement(self):
        basis_tensor = np.kron(self.basis_state_one,self.basis_state_two)
        H_tensor_I = np.kron(self.hadamard_gate,self.I)
        CNOT_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        a1 = np.dot(H_tensor_I,basis_tensor)
        bell_state = (self.hadamard_constant * 2)*np.dot(CNOT_gate,a1)
        return bell_state

    #write functions to clean up reused code. Looks sloppy and the same

    def equation_builder(self,circuit_layout):
        equation = []
        for i in range(len(circuit_layout[0])):
            if i == 0:
                basis_tensor = ()
                if circuit_layout[0][i]+circuit_layout[0][i+1] == 'e0':
                    if circuit_layout[1][i]+circuit_layout[1][i+1] == 'e0':
                        basis_tensor = (np.array([[1],[0]]),np.array([[1],[0]]))
                    elif circuit_layout[1][i]+circuit_layout[1][i+1] == 'e1':
                        basis_tensor = (np.array([[1],[0]]),np.array([[0],[1]]))
                elif circuit_layout[0][i]+circuit_layout[0][i+1] == 'e1':
                    if circuit_layout[1][i]+circuit_layout[1][i+1] == 'e0':
                        basis_tensor = (np.array([[0],[1]]),np.array([[1],[0]]))
                    elif circuit_layout[1][i]+circuit_layout[1][i+1] == 'e1':
                        basis_tensor = (np.array([[0],[1]]),np.array([[0],[1]]))                    
                equation.append(basis_tensor)
            if circuit_layout[0][i] == '[':
                tensor = ()
                if circuit_layout[0][i+1] == 'I':
                    if circuit_layout[1][i+1] == 'I':
                        tensor = (self.I,self.I)
                    elif circuit_layout[1][i+1] == 'H':
                        tensor = (self.I,self.hadamard_gate)
                    elif circuit_layout[1][i+1] == 'X':
                        tensor = (self.I,self.NOT_gate)
                elif circuit_layout[0][i+1] == 'H':
                    if circuit_layout[1][i+1] == 'I':
                        tensor = (self.hadamard_gate,self.I)
                    elif circuit_layout[1][i+1] == 'H':
                        tensor = (self.hadamard_gate,self.hadamard_gate)
                    elif circuit_layout[1][i+1] == 'X':
                        tensor = (self.hadamard_gate,self.NOT_gate)
                elif circuit_layout[0][i+1] == 'X':
                    if circuit_layout[1][i+1] == 'I':
                        tensor = (self.NOT_gate,self.I)
                    elif circuit_layout[1][i+1] == 'H':
                        tensor = (self.NOT_gate,self.hadamard_gate)
                    elif circuit_layout[1][i+1] == 'X':
                        tensor = (self.NOT_gate,self.NOT_gate)  
                equation.append(tensor)
            if circuit_layout[0][i] == '|':
                if circuit_layout[0][i+1] == 'C':
                    equation.append(self.CNOT_gate)
                if circuit_layout[0][i+1] == 'A':
                    equation.append(self.AT_gate)
                if circuit_layout[0][i+1] == 'S':
                    equation.append(self.SWAP_gate)
        return equation


    def evaluate_equation(self,circuit_layout):
        equation = self.equation_builder(circuit_layout)
        new_equation = []
        output = 0
        for i in range(len(equation)):
            if type(equation[i]) is tuple:
                new_equation.append(np.kron(equation[i][0],equation[i][1]))
            else:
                new_equation.append(equation[i])

        #for i in range(len(new_equation)):
            #output = np.dot()
        return new_equation




test_string = ["e0---[I]---[H]---|CNOT|---[H]---(e0)",
               "e0---[X]---[H]---|CNOT|---[H]---z"]


e0 = np.array([[1],[0]])
e1 = np.array([[0],[1]])
CNOT_gate = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
SWAP_gate = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
AT_gate = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
myCircuit = QCircuit(e0,e0,CNOT_gate)

output = myCircuit.deutsch_algorithm()

e10 = np.kron(e1,e0)
e11 = np.kron(e1,e1)
hadamard_constant = 1/(2**(1/2))
superposition = hadamard_constant*(e10-e11)

print("Quantum circuit layout for Deutsch algorithm:")
print()
print("e0---[I]---[H]---|  |---[H]---(e0)")
print("                 |Uf|")
print("e0---[X]---[H]---|  |---[I]---z")
print()
print("Input boolean function gate Uf: ")
print(CNOT_gate)
print("Quantum circuit evaluated with numeric linear algebra: ")
print(output)
print("Quantum circuit evaluated with properties of tensor product spaces: ")
print(superposition)


bell_state = myCircuit.qbit_entanglement()

print()
print("Quantum circuit layout for 2 qbit entanglement:")
print()
print("e0---[H]--*---")
print("          |")
print("e0---[I]--@---")
print()
print("Output state (Bell state): ")
print(bell_state)


print()

print(myCircuit.evaluate_equation(["e0---[I]---[H]---|CNOT|---[H]---(e0)","e0---[X]---[H]---|CNOT|---[I]---z"]))



print()



def bin_num(num):
    return int(num[1]) + 2*int(num[0])

def CV():
    Im = [[1,0],[0,1]]
    V = [[1+1j,1-1j],[1-1j,1+1j]]
    matrix = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

    permutations = ['000','001','010','011','100','101','110','111']
    for i in permutations:
        a = i[0]
        b = i[1]
        c = i[2]
        matrix[bin_num('0' + a)][bin_num(b+c)] = Im[int(b)][0] * Im[int(a)][int(c)]
        matrix[bin_num('1' + a)][bin_num(b+c)] = Im[int(b)][1] * V[int(a)][int(c)]

    return matrix

print(np.array(CV()))

def bin_num_3(num):
    return 4*int(num[0]) + 2*int(num[1]) + int(num[2])

def CIV():
    Im = [[1,0],[0,1]]
    V = [[1+1j,1-1j],[1-1j,1+1j]]
    matrix = [[0 for i in range(8)] for j in range(8)]

    permutations = ['{0:0{width}b}'.format(v, width=5) for v in range(2**5)]
    for i in permutations:
        a = i[0]
        b = i[1]
        c = i[2]
        d = i[3]
        e = i[4]
        matrix[bin_num_3('0' + a + b)][bin_num_3(c + d + e)] = Im[0][int(c)]*Im[int(a)][int(d)]*Im[int(b)][int(e)]
        matrix[bin_num_3('1' + a + b)][bin_num_3(c + d + e)] = Im[1][int(c)] * Im[int(a)][int(d)] * V[int(b)][int(e)]

    return matrix

print()

CIV = CIV()
for i in range(len(CIV)):
    if CIV[i][i] == 1:
        CIV[i][i] = 2
CIV = np.array(CIV)

print("1/2*"+str(np.array(CIV)))
print()

CIV = np.array(CIV)
I = np.array([[1,0],[0,1]])
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CV = np.array([[2,0,0,0],[0,2,0,0],[0,0,1+1j,1-1j],[0,0,1-1j,1+1j]])
CVdagger = np.transpose(np.array([[2,0,0,0],[0,2,0,0],[0,0,1-1j,1+1j],[0,0,1+1j,1-1j]]))

ItCV = np.kron(I,CV)
CNOTtI = np.kron(CNOT,I)
ItCVdagger = np.kron(I,CVdagger)

a1 = np.dot(CNOTtI,CIV) #1/2 from CIV
a2 = np.dot(ItCVdagger,a1) #1/2 from CVdagger
a3 = np.dot(CNOTtI,a2)
a4 = np.dot(ItCV,a3) #1/2 from CV

print((1/2)*(1/2)*(1/2)*a4)

print()

l1 = np.kron(e0,e0)
l3 = np.dot(CNOT,CNOT)
l4 = np.dot(l3,CNOT)
l5 = np.kron(l4,e0)

e0 = np.array([[1],[0]])
e0dagger = np.array([[1,0]])
e0e0dagger = np.dot(e0,e0dagger)
l6 = np.kron(CNOT,e0e0dagger)
l7 = np.dot(l6,np.dot(l6,l6))
print(l6)
    
H = np.array([[1,1],[1,-1]])
print(CV * CV)
print(CV)

CT = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1+1j]])
print(CT*CT)

print(np.kron(I,H))

ItH = np.kron(I,H)
CS = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1j]])

a1 = np.dot(CS,ItH)
a2 = np.dot(a1,CS)
a3 = np.dot(a2,ItH)

print(a3)

a = np.array([[1,2,4,2,1]])
b = np.transpose(a)

print(np.dot(b,a))

a = np.array([[1,0],[0,-1]])
print(a**2)

b = np.array([[0,0],[1,2]])

print(np.dot(b,a))
     