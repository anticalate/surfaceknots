import sys, time, json, pickle
from surfaceknots.fast_enum import enumerate_targeted_prime_ch_census
from surfaceknots.moves import has_certified_cyclic_unknotted_tube
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
ch=int(sys.argv[1]); nw=int(sys.argv[2])
t=time.time()
c,st=enumerate_targeted_prime_ch_census(ch=ch,m_min=1,m_max=ch-1,target_components=0,target_orientability=-1,num_workers=nw,filter_cyclic_tubes=False)
print('ch',ch,'kept',len(c),st,round(time.time()-t,1),flush=True)
out=[]
for d in c:
    topo=d.surface_topology()
    flag=has_certified_cyclic_unknotted_tube(d)
    wp=wirtinger_presentation(d)
    fq=finite_quotient_profile(wp,("S3","Q8","A4","S4","A5"))
    rec=dict(epd=d.to_epd(),ncomp=topo.num_components,types=[cc.homeomorphism_type for cc in topo.components],e=[cc.normal_euler for cc in topo.components],tube=flag,hom=dict(fq),ngen=wp.num_generators)
    out.append(rec)
    print(('TUBE ' if flag else '     '), rec['types'], rec['e'], rec['hom'], 'gens',wp.num_generators, flush=True)
json.dump(out,open(f'{sys.argv[3]}/tube_ch{ch}.json','w'),indent=1)
