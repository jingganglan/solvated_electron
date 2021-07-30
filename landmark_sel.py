import numpy as np
idx_landmark=np.genfromtxt('./cursel.landmarks')
idx_landmark=sorted(idx_landmark)

def get_frame_size(input_name):
    counter=0
    with open(input_name, 'r') as infile:
        for line in infile:
            counter=counter+1
            if "begin" in line:
                idx_begin = counter
            if "end" in line:
                idx_end = counter
                break
    return idx_end - idx_begin + 1



def read_frame(input_name,idx_frame):
    counter = 0
    frame_size = get_frame_size(input_name)
    frame_content = ""
    
    with open (input_name, 'r') as infile:
        
        for line in infile:

            if ( counter >= frame_size*idx_frame ):
                
                frame_content += line
            
            if ( counter == frame_size*(idx_frame+1)-1):
                
                break
                
            counter = counter + 1
            
    return frame_content

def frame2xyz(frame_content,i):
    lines = frame_content.splitlines()
    n_atom = len(lines) - 8
    cell_a=lines[2:5][0].split()[1:4]
    cell_b=lines[2:5][1].split()[1:4]
    cell_c=lines[2:5][2].split()[1:4]
    f = open("mp2_"+str(int(i))+".xyz", "w")
    f.write(str(n_atom)+"\n")
    f.write(("ABC"+"  "+str(cell_a[0])+"  "+str(cell_b[1])+"  "+str(cell_c[2])+"  "+"#"+str(i))+"  "+"\n")
    for i in range(n_atom):
        atom  = lines[5:5+n_atom][i].split()[4]
        coord = np.array(lines[5:5+n_atom][i].split()[1:4],dtype='float')
        f.write(str(atom)+"  "+str(coord[0]*0.529177)+"  "+str(coord[1]*0.529177)+"  "+str(coord[2]*0.529177)+"\n")
    
    f.close()

    
for i in idx_landmark:
    frame_content = read_frame("./input.data",i)
    frame2xyz(frame_content,i)

