#!/usr/bin/env python
from subprocess import Popen, PIPE

conf_level = 0.95
upper = 1e-3
lower = -1e-3

val_cache = {}

def find_limit(value,clevel):
    alpha = 1-clevel
    if value not in val_cache:
        combine_args = ['combine',
                        'Example_AQGC2_workspace.root',
                        '-M','HybridNew',
                        '--freq',
                        '--testStat','PL',
                        '--rule','CLsplusb',
                        '--toysH','500',
                        '--clsAcc','0.002',                        
                        '--singlePoint','a0W=0,aCW=%e'%value]
        
        result = Popen(combine_args,stdout=PIPE,stderr=PIPE).communicate()[0]
        result = result.split('\n')[-3]
        result = result.split('=')[-1]
        result = result.split('+/-')
        central_val = float(result[0])
        error       = float(result[1])
        val_cache[value] = [central_val,error]    
       
    print 'limit at %e: %f +/- %e'%(value,
                                    val_cache[value][0],
                                    val_cache[value][1])
    return val_cache[value][0] >= alpha, val_cache[value]
    
def get_intervals(lower,guess,upper,conf_level):
    # check midpoint    
    up_inc, up_and_err = find_limit(upper,conf_level)
    lo_inc, lo_and_err = find_limit(lower,conf_level)
    gs_inc, gs_and_err = find_limit(guess,conf_level)    
    avg_err = (up_and_err[1] + lo_and_err[1] + gs_and_err[1])/ 3.0
    diff_lo_gs = abs(lo_and_err[0] - gs_and_err[0])
    diff_up_gs = abs(up_and_err[0] - gs_and_err[0])
    print diff_lo_gs, diff_up_gs, avg_err
    if ( diff_lo_gs < avg_err and diff_up_gs < avg_err and
         lo_inc and gs_inc and up_inc ):
        return guess
    else:
        # neither of regions left of upper contain the limit
        # reset midpoint to be between guess and upper
        if not lo_inc and not gs_inc and up_inc:
            midpoint = 0.5*(upper+guess)
            return get_intervals(guess,midpoint,upper,conf_level)
        # upper and guess are both in the limit
        # search between lower and guess for the limit
        if not lo_inc and gs_inc and up_inc:
            midpoint = 0.5*(lower+guess)
            return get_intervals(lower,midpoint,guess,conf_level)
        # neither of the regions right of lower contain the limit
        # search for the limit between lower and guess
        if lo_inc and not gs_inc and not up_inc:
            midpoint = 0.5*(lower+guess)
            return get_intervals(lower,midpoint,guess,conf_level)
        # both lower and guess are in the limit
        # search for the limit between guess and upper
        if lo_inc and gs_inc and not up_inc:
            midpoint = 0.5*(upper+guess)
            return get_intervals(guess,midpoint,upper,conf_level)
        # there are no points which are inside the limits
        # fork into two!
        if not lo_inc and not gs_inc and not up_inc:
            midpoint_up = 0.5*(upper+guess)
            midpoint_lo = 0.5*(lower+guess)
            val_up = get_intervals(guess,midpoint_up,upper,conf_level)
            val_lo = get_intervals(lower,midpoint_lo,guess,conf_level)
            return [val_lo,val_up]
        # there are the middle is in and the endpoints are out
        # fork into two!
        if not lo_inc and gs_inc and not up_inc:
            midpoint_up = 0.5*(upper+guess)
            midpoint_lo = 0.5*(lower+guess)
            val_up = get_intervals(guess,midpoint_up,upper,conf_level)
            val_lo = get_intervals(lower,midpoint_lo,guess,conf_level)
            return [val_lo,val_up]
        # there are the middle is in and the endpoints are out
        # fork into two!
        if lo_inc and not gs_inc and up_inc:
            midpoint_up = 0.5*(upper+guess)
            midpoint_lo = 0.5*(lower+guess)
            val_up = get_intervals(guess,midpoint_up,upper,conf_level)
            val_lo = get_intervals(lower,midpoint_lo,guess,conf_level)
            return [val_lo,val_up]
        
   
if __name__ == '__main__':
    included = find_limit(6.2159e-4,conf_level)
    print included
    
    included = find_limit(0.000625610351562,conf_level)
    print included
    
    included = find_limit(-0.000625610351562,conf_level)
    print included

    result = get_intervals(lower,0.5*(upper+lower),upper,conf_level)
    print result

