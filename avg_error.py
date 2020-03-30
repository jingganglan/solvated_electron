import numpy as np

def confidence(inf,size,conf):
    n_element = int(inf.shape[0]/size)
    avg=np.zeros(n_element)
    top=np.zeros(n_element)
    bot=np.zeros(n_element)
    run_avg=np.zeros(n_element)
    trust = 0.5
    a=np.std(inf)
    while trust < conf:
        for n in range(n_element):
            avg[n]=np.average(inf[n*size:(n+1)*size])
            run_avg[n]=np.average(avg[:n])
            top[n]=run_avg[n]+a/np.sqrt(n)
            bot[n]=run_avg[n]-a/np.sqrt(n)

        list_a=np.where( (avg < top ) & (avg > bot) )
        trust=np.size(list_a)/inf.shape[0]*size
        
        a = a+a*(conf - trust)
    print(trust,a,"average",run_avg[n], "error:", a/np.sqrt(n_element))

    return avg, top, bot

