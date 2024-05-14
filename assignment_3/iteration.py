#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

def iterations(x_range, y_range, step):
    """
    A function 
    """
    abs_lst = []
    div_lst = []
    iter_lst = []
    for x in np.arange(x_range[0], x_range[1], step):
        for y in np.arange(y_range[0], y_range[1], step):
            c = x + y * 1j
            z = 0
            iter = 0
            
            for _ in range(10):
                z = z**2 + c
                iter += 1
                
                if abs(z) > 2:
                    div_lst.append(c)
                    iter_lst.append(iter)
                    break
            else:
                abs_lst.append(c)

    return abs_lst, div_lst, iter_lst

