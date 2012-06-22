'''
test anova with Data object
'''

import numpy as np

from pyCMBS import *

from pylab import *

from anova import *




a = ANOVA()

#~ a.add_experiment('e2')


#~ d2=Data(None,None)
#~ d3=Data(None,None)
#~ d4=Data(None,None)


#~ x=rand(10,100,200)
#~ d1.data = np.ma.array(x,mask=x>0.9)
#~ a.add_data('e2',d1)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d2.data = np.ma.array(x,mask=x>0.9)
#~ a.add_data('e2',d2)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d3.data = np.ma.array(x,mask=x>2.)
#~ a.add_data('e1',d3)
#~ del x
#~
#~ x=rand(10,100,200)
#~ d4.data = np.ma.array(x,mask=x>2.)
#~ a.add_data('e1',d4)


#http://people.richland.edu/james/lecture/m170/ch13-2wy.html
a.add_experiment('e1')
a.add_experiment('e2')
a.add_experiment('e3')

#ens #1
d1=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [106,95,94,103,100]
d1.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e1',d1)

d2=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [110,98,100,108,105]
d2.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e2',d2)

d3=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [94,86,98,99,94]
d3.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e3',d3)

#ens #2
d4=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [110,100,107,104,102]
d4.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e1',d4)

d5=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [112,99,101,112,107]
d5.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e2',d5)

d5=Data(None,None)
x = np.zeros((5,1,1))
x[:,0,0] = [97,87,99,101,98]
d5.data = np.ma.array(x,mask=np.isnan(x))
a.add_data('e3',d5)



#~ stop
#~ groups = [
          #~ [[106,95,94,103,100],
          #~ [110,98,100,108,105],
          #~ [94,86,98,99,94]],
#~
          #~ [[110,100,107,104,102],
          #~ [112,99,101,112,107],
          #~ [97,87,99,101,98]]
#~
          #~ ]




a.analysis()
