from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/2W83.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta






# SAS points

load "data/2W83.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [4555,3888,3740,3873,4554,3887,3732,3890,3729,3884,3878,3736,3770,3746,3772,3747,3771,3737,3739,3741,3738,4068,3764,3894,3763,3895,3748,3750,3758] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.310,0.702]
select surf_pocket2, protein and id [103,112,100,117,99,93,89,91,90,92,226,124,225,116,243,237,240,111,421,82,85] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.498,0.361,0.902]
select surf_pocket3, protein and id [1527,1528,1414,1418,1419,1426,2490,1549,1720,1545,1395,1532,1542,1533,1539,2210,1392,1394,2209,1391,1393,1550,1405,1402,1413,1401,1723,1387,1732,1727] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.525,0.278,0.702]
select surf_pocket4, protein and id [1177,1188,925,1187,924,898,899,100,905,908,907,906,92,125,226,124] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.851,0.361,0.902]
select surf_pocket5, protein and id [4835,4555,4553,4554,4552,4548,3770,3772,3747,3771,4546,3739,4834,4571] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.600]
select surf_pocket6, protein and id [503,502,790,504,811,809,3374,810,517,501,515,3366,3342,2832,3341,3340,3311,3313,2803,2802,2801] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.596]
select surf_pocket7, protein and id [3473,294,3475,273,280,3450,276,295,283,2943,2967,278,2911,2934,2968] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.325]
select surf_pocket8, protein and id [3466,3436,3465,3412,2977,1597,1581,1454,1452,3444,1585] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.478,0.361]
select surf_pocket9, protein and id [1557,1559,1758,1759,1560,3382,1756,1562,1784] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.510,0.278]
select surf_pocket10, protein and id [4615,4563,4631,968,986,984,4624,976,977,967] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.831,0.361]
select surf_pocket11, protein and id [4803,4882,4883,985,989,992,993,4905,4881,4907,4787,4912,4910,1030,991,1022] 
set surface_color,  pcol11, surf_pocket11 




deselect

orient
