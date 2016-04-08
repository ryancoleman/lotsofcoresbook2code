name='diamond333_hch'

calc = Calculator(name +'.gpw', 
         kpts=(6,6,6),txt=name + '_rec.txt')
calc.set_positions()

r = RecursionMethod(calc)
r.run(600) 

r.run(1400, 
      inverse_overlap='approximate')

r.write(name + '_600_1400a.rec',
        mode='all')
