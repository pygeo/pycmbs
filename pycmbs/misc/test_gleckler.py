from pyCMBS import *
G = GlecklerPlot()
#register first models
G.add_model('mod1')
G.add_model('mod2')
G.add_model('mod3')
G.add_model('mod4')
G.add_model('mod5')
#then register variables
G.add_variable('ta')
G.add_variable('P')
#after that you can add values to be plotted; pos=1 mean that result is plotted in upper triangle
G.add_data('ta','mod1',1.,pos=1)
G.add_data('ta','mod2',2.,pos=1)
G.add_data('ta','mod3',3.,pos=1)
G.add_data('ta','mod4',4.,pos=1)
G.add_data('ta','mod5',5.,pos=1)

#another dataset
G.add_data('ta','mod1',0.5,pos=2)
G.add_data('ta','mod2',0.4,pos=2)
G.add_data('ta','mod3',0.3,pos=2)
G.add_data('ta','mod4',0.2,pos=2)
G.add_data('ta','mod5',0.1,pos=2)


G.add_data('P','mod1',0.25,pos=1)
G.add_data('P','mod2',-0.25,pos=2)
G.add_data('P','mod3',-0.25,pos=1)

G.write_ranking_table('ta', 'test_rnk_ta.tex')
G.write_ranking_table('ta', 'test_rnk_ta.md', fmt='markdown')

G.write_ranking_table('P', 'test_rnk_P.tex')

G.plot(show_value=True) #do plot
