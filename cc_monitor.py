import sys, re, os, colorsys, random, platform
from utilities import get_list_from_file, write_list_to_file
if platform.system()=='Windows':
    import cairo
    workdir = 'C:\\Users\\miken\\Dropbox\\Grad School\\Phylogenetics\\work\\test-queue\\'
    file_list = os.path.join(workdir,'text_files_completed.txt')
else:
    workdir = '/projects/tallis/nute/qdumps'

import math

def get_date_from_file(filename):
    fo, fi = os.path.split(filename)
    a = fi.replace('.txt','')
    a = a.replace('CDT_','')
    fn_date_str = a
    b= a.split('_')
    # print b
    # a = a.replace('_',' ')
    date_str =' '.join(b[2:4]) + ' ' + b[7] +' -- ' + ':'.join(b[4:7])
    return date_str, fn_date_str

def get_folder_inventory():
    text_files_done=get_list_from_file(file_list)
    text_files_present = os.listdir(os.path.join(workdir,'text'))
    new_files = list(set(text_files_present).difference(text_files_done))
    return new_files, text_files_done


def process_file(filename=None):
    if filename==None:
        filename = os.path.join(workdir,'text','queue_jobs_Sep_9_09_53_35_CDT_2016.txt')
    else:
        filename = os.path.join(workdir,'text',filename)

    myf = open(filename,'r')
    prev = myf.readline()
    jobs={}
    nodes={}
    doing_jobs = True
    idre = re.compile("^Job Id: (?P<id>[0-9]+).*")

    myargs = {}
    mainkey = ''
    lastkey = ''

    jct = 0
    nct = 0

    for line in myf:
        if len(line.strip())==0:
            if mainkey<>'' and len(myargs)>0:
                if doing_jobs==True:
                    jobs[mainkey]=myargs
                    # jct += 1
                    # print "job %s: %s" % (jct,mainkey)
                else:
                    nodes[mainkey]=myargs
                    # nct += 1
                    # print "node %s: %s" % (nct,mainkey)
                mainkey = ''
                lastkey = ''
                myargs = {}

        elif line[0] == '\t':
            if lastkey<>'':
                myargs[lastkey] += line[1:].strip()
        elif line[0:4] == '    ':
            a = line[4:].split(' = ')
            myargs[a[0].strip()] = a[1].strip()
            lastkey = a[0]
        else:
            temp_key = line.strip()
            if temp_key[0:3]<>'Job':
                if temp_key[0:3] in ['tau','gol']:
                    doing_jobs = False
                    # print temp_key

            if doing_jobs == True:
                try:
                    mainkey = idre.match(temp_key).group('id')
                except:
                    print line
                #     print 'error...'
                #     return None
            else:
                mainkey = line.strip()

    ds, ds_fi = get_date_from_file(filename)
    return nodes, jobs, ds, ds_fi


def color_scale_set(total):
    lower_lim = .40

    side = int(float(total) ** .333 + 1)
    gap = 1.0 / float(side)

    ineliglbe_pts = int((lower_lim - gap * 5 / 6) / gap + 1) ** 3
    ineliglbe_pts += 1 #white is not allowed

    while side ** 3 - ineliglbe_pts < total:
        side += 1
        gap = 1.0 / float(side)
        ineliglbe_pts = int((lower_lim - gap * 5 / 6) / gap + 1) ** 3

    coords = []
    for i in range(side):
        c = gap * 5 / 6 + gap * float(i)
        coords.append(c)

    locus = []
    backrange = range(side)
    backrange.sort(reverse=True)
    for i in backrange:
        for j in backrange:
            for k in backrange:
                newc = (coords[i], coords[j], coords[k])
                if max(newc) > lower_lim and min(i,j,k)<max(backrange):
                    locus.append(newc)

    sss = side * side * side - ineliglbe_pts
    perm = get_ideal_permutation(sss)
    finals = []
    for i in perm:
        finals.append(locus[len(perm) - i - 1])
    return finals


def get_ideal_permutation(els):
    random.seed(100)
    m1 = int(els / 2)
    m2 = m1 + 1
    Ai = range(1, m1)
    Bi = range(m2 + 1, els + 1)
    A = random.sample(Ai, len(Ai))
    B = random.sample(Bi, len(Bi))
    out = []
    out.append(m1 - 1)
    if len(B) <> len(A):
        out.append(B.pop() - 1)
    for i in range(len(A)):
        out.append(A.pop() - 1)
        out.append(B.pop() - 1)
    out.append(m2 - 1)
    return out

def make_graphic(nodes,jobs,date_str, fi_date_str,graphicpath=None):
    jre = re.compile('([0-9]{7})')
    if graphicpath==None:
        graphic_name_spec = 'graphic_' + fi_date_str+ '.png'
        graphicpath_spec=os.path.join(workdir,'visuals',graphic_name_spec)
        graphicpath_gen = os.path.join(workdir, 'graphic.png')
    mycolors = color_scale_set(100)


    WIDTH, HEIGHT = 1600, 1200
    cell_wd, cell_ht = 60, 18
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0,0,WIDTH,HEIGHT)
    ctx.fill()
    userlist=[]
    ctx.translate(20, 20)

    ctx.move_to(20,20)
    ctx.set_font_size(12)
    ctx.set_source_rgb(0, 0, 0)
    # ctx.fill()
    index = 0
    for i in nodes.keys():
        if 'jobs' in nodes[i].keys():
            joblist = jre.findall(nodes[i]['jobs'])
            if len(joblist)>1:
                print i + ' has %s jobs' % len(joblist)
            if 'euser' in jobs[joblist[0]].keys():
                nodes[i]['occupant'] = jobs[joblist[0]]['euser']
                userlist.append(nodes[i]['occupant'])
            nodes[i]['joblist'] = joblist

    userset = list(set(userlist))
    users_on_nodes = {}
    for i in userset:
        users_on_nodes[i]=[]
    for i in nodes.keys():
        if 'occupant' in nodes[i].keys():
            users_on_nodes[nodes[i]['occupant']].append(i)
        else:
            s=nodes[i]['state']
            if s not in users_on_nodes.keys():
                users_on_nodes[s]=[]
            users_on_nodes[s].append(i)

    temp = {}
    for i in userset:
        a=float(len(users_on_nodes[i]))
        if a in temp.keys():
            b=a+random.random()
        else:
            b=a+0.0
        temp[b]=i
    counts = list(temp.keys())
    counts.sort(reverse=True)

    index = 0
    for i in counts:
        usr = temp[i]
        usr_nodes = users_on_nodes[usr]
        for j in usr_nodes:
            row = float(index % 50)
            col = float(int(index/50))
            clr = mycolors[userset.index(usr)]
            ctx.set_source_rgb(*clr)
            ctx.rectangle((col) * cell_wd, (row - 1) * cell_ht + 5, cell_wd, cell_ht)
            ctx.fill()
            ctx.set_source_rgb(0, 0, 0)
            ctx.move_to(col*cell_wd,row*cell_ht)
            ctx.show_text(usr)
            index+=1

    otherstates=set(users_on_nodes.keys()).difference(set(userlist))
    if len(otherstates)>0:
        for i in otherstates:
            if i=='free':
                clr = (1,1,1)
            else:
                clr = (0,0,0)
            for j in users_on_nodes[i]:
                row = float(index % 50)
                col = float(int(index / 50))
                ctx.set_source_rgb(*clr)
                ctx.rectangle((col) * cell_wd, (row - 1) * cell_ht + 5, cell_wd, cell_ht)
                ctx.fill()
                ctx.set_source_rgb(0, 0, 0)
                ctx.move_to(col * cell_wd, row * cell_ht)
                ctx.show_text(j)
                index += 1

    ypos=20
    vertical_border = 5
    ctx.move_to((col+1)*cell_wd+60,ypos)
    ctx.set_font_size(20)
    ctx.show_text('Jobs Running by User:')
    (x, y, wd, ht, dx, dy) = ctx.text_extents('Jobs Running by User:')
    ctx.set_font_size(12)
    left_edge = (col+1)*cell_wd+35+75
    # right_edge = WIDTH-100
    right_edge = left_edge + 300
    ctx.move_to(right_edge,ypos)
    ctx.set_source_rgb(.25, .25, .25)
    ctx.show_text('100')
    ctx.fill()
    ypos += 10
    scale = (float(abs(right_edge-left_edge))/100)
    # ctx.set_source_rgb(.25,.25,.25)


    index = 0
    user_list_order=[]
    for i in counts:
        usr = temp[i]
        user_list_order.append(usr)
        ct = len(users_on_nodes[usr])
        (x, y, wd, ht, dx, dy) = ctx.text_extents(usr)
        ypos+=ht + vertical_border
        ctx.move_to(left_edge - wd - 5,ypos)
        ctx.show_text(usr)
        ctx.move_to(left_edge,ypos)
        ctx.rectangle(left_edge,ypos -ht ,scale*ct,ht)
        ctx.fill()

    nouser_count = 0
    user_queued_job_counts={}
    user_queued_node_request_counts={}
    users_queued_only=[]

    for i in jobs.keys():
        if 'job_state' in jobs[i].keys() and jobs[i]['job_state']=='Q':
            if 'queue' in jobs[i].keys() and jobs[i]['queue'] == 'secondary':
                if 'euser' in jobs[i].keys():
                    usr = jobs[i]['euser']
                    if usr not in user_queued_job_counts.keys():
                        user_queued_job_counts[usr]=0
                    user_queued_job_counts[usr]+=1
                else:
                    nouser_count+=1
                if 'job_array_request' in jobs[i].keys():
                    arrstr = jobs[i]['job_array_request'].split('-')
                    nd_r_ct = int(arrstr[1])
                else:
                    nd_r_ct = 1

                if usr not in user_queued_node_request_counts.keys():
                    user_queued_node_request_counts[usr] = 0
                user_queued_node_request_counts[usr] += nd_r_ct
    users_running_only = set(userlist).difference(set(user_queued_node_request_counts.keys()))
    users_queued_only = list(set(user_queued_node_request_counts.keys()).difference(set(userlist)))
    for i in users_running_only:
        user_queued_node_request_counts[i]=0

    queed_only_cts = {}
    for i in users_queued_only:
        ct = float(user_queued_node_request_counts[i])
        if ct in queed_only_cts.keys():
            ct=ct + random.random()
        queed_only_cts[ct]=i
    queed_only_cts_list = list(queed_only_cts.keys())
    queed_only_cts_list.sort(reverse=True)


    ypos = 20 + 16*(ht + vertical_border)
    ctx.move_to(left_edge + 100, ypos)
    ctx.set_font_size(20)
    # (x, y, wd, ht, dx, dy) = ctx.text_extents('# axci')
    # print 'x: %s, y: %s, wd: %s, ht: %s, dx: %s, dy: %s' % (x, y, wd, ht, dx, dy)
    ctx.show_text('# Jobs Queued by User:')
    (x, y, wd, ht, dx, dy) = ctx.text_extents('# Jobs Queued by User:')
    # print 'x: %s, y: %s, wd: %s, ht: %s, dx: %s, dy: %s' % (x, y, wd, ht, dx, dy)
    ypos += 10
    ctx.set_font_size(12)
    left_edge = left_edge + 100
    # right_edge = WIDTH-100
    right_edge = left_edge + 300
    scale = (float(abs(right_edge - left_edge)) / 100)
    ctx.set_source_rgb(.25, .25, .25)
    ctx.move_to(right_edge,ypos-12)
    ctx.show_text('100')
    ctx.fill()

    total_order = user_list_order + users_queued_only
    for i in total_order:
        ct = user_queued_node_request_counts[i]
        (x, y, wd, ht, dx, dy) = ctx.text_extents(i)
        ypos+=10 + vertical_border
        ctx.move_to(left_edge - wd - 5,ypos)
        ctx.show_text(i)
        ctx.move_to(left_edge,ypos)
        ctx.rectangle(left_edge,ypos -10+1 ,scale*ct,9)
        ctx.fill()

    ctx.move_to(25,50*cell_ht+35)
    ctx.set_font_size(20)
    ctx.show_text('As of:  ' + date_str)
    ctx.fill()

    surface.write_to_png(graphicpath_spec)  # Output to PNG
    surface.write_to_png(graphicpath_gen)  # Output to PNG


if __name__=='__main__':
    files_to_do, already_done = get_folder_inventory()
    for i in files_to_do:
        print 'running %s' % i
        n, j, dstr, f_dtstr = process_file(i)
        make_graphic(n,j, dstr, f_dtstr)
        already_done.append(i)
    write_list_to_file(already_done,file_list)


    pass