#!/usr/bin/python
def list_of_groups(init_list, childern_list_len):
  list_of_groups = zip(*(iter(init_list),) *childern_list_len)
  end_list = [list(i) for i in list_of_groups]
  count = len(init_list) % childern_list_len
  end_list.append(init_list[-count:]) if count !=0 else end_list
  return end_list
import re
import os
#########################################################################
total_list = []
for filename in os.listdir('.'):
#for filename in ['common.h']:
  if not ('.c' in filename or '.h' in filename):
    continue;
  print 'processing ' + filename
  var_param = r"(struct)?\s*([\w\[0-9\]]+)\s*([*&]*)\s*([\w\[0-9\]]*)"
  #pattern = r"(\w+)\s+[*,&]?\s?(\w+)\s?+\("
  pattern1 = r"\s*%s\s+%s\s*\(\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param)
  pattern2 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param)
  pattern3 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param)
  pattern4 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param)
  pattern5 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern6 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern7 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern8 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern9 = r"\s*%s\s+%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s*([;\{])" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pat1 = re.compile(pattern1)
  pat2 = re.compile(pattern2)
  pat3 = re.compile(pattern3)
  pat4 = re.compile(pattern4)
  pat5 = re.compile(pattern5)
  pat6 = re.compile(pattern6)
  pat7 = re.compile(pattern7)
  pat8 = re.compile(pattern8)
  pat9 = re.compile(pattern9)
  pattern = [pat1, pat2, pat3, pat4, pat5, pat6, pat7, pat8, pat9]
  f = open(filename, "r")
  line_index = 1
  space_pat = re.compile(r'}\s*')
  segment = []
  begin = []
  end = []
  flag = False
  rep_info = []
  head = False
  new_content = ''
  for line in f.readlines():
    print '\n\n'
    print line
    for (index, pat) in zip(range(0,9), pattern):
      res = pat.match(line)
      if res:
        list_res = list_of_groups(res.groups(), 4)
        begin = [line_index, index, list_res]
        print begin
        rep_info = []
        for (ind_t, t) in zip(range(0,len(list_res[2:])), list_res[2:]):
          if t[0] == 'struct' and t[2] == '' and t[3] == '':
            line = line.replace(' ' + t[1] + ',', ' ' + t[1] + '*,');
            line = line.replace(' ' +  t[1] + ')', ' ' + t[1] + '*)');
            rep_info.append((1, t[1]))
            #total_list.append([list_res[1][1], ind_t]) 
          if t[0] == 'struct' and t[2] == '' and t[3] != '':
            rep_info.append((3, t[3]))
            line = line.replace(' ' + t[3] + ',', ' *' + t[3] + ',');
            line = line.replace(' ' + t[3] + ')', ' *' + t[3] + ')');
            if list_res[-1][0] == '{':
              total_list.append([list_res[1][1], ind_t, filename, line_index, 0, t[3]]) 
              flag = True
              head = True
        print rep_info    
    res2 = space_pat.match(line)
    if res2:
      end = [line_index, res2.groups()]
      print  end
      segment.append([begin, end])
      print flag
      if flag:
        for (ind_tt, res) in zip(range(0, len(total_list)), total_list):
          if res[0] == total_list[-1][0] and res[2] == total_list[-1][2]:
            total_list[ind_tt][4] = line_index
      flag = False
    if flag:
      if head:
        head = False
      else:
        for t in rep_info:
          if not t == '':
            line = line.replace(t[1] + '.', t[1] + '->');
    new_content += line
    line_index += 1
  print segment
  f.close()
  f2 = open(filename, 'w')
  f2.write(new_content)
  f2.close()
new_list = []
for i in total_list:
  if not i in new_list:
    new_list.append(i)
print new_list

#########################################################################
#for filename in ['macdrp.c']:
for filename in os.listdir('.'):
  if not ('.c' in filename or '.h' in filename):
    continue;
  print 'processing ' + filename
  var_param = r"([*&]*\s*[\w\[0-9\]\+\-\*\/]+)"
  #pattern = r"(\w+)\s+[*,&]?\s?(\w+)\s?+\("
  pattern1 = r"(\s*)%s\s*\(\s*%s\s*\)\s?;" % (var_param, var_param)
  pattern2 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param)
  pattern3 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param)
  pattern4 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param)
  pattern5 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param, var_param)
  pattern6 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern7 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern8 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pattern9 = r"(\s*)%s\s*\(\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*,\s*%s\s*\)\s?;" % (var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param, var_param)
  pat1 = re.compile(pattern1)
  pat2 = re.compile(pattern2)
  pat3 = re.compile(pattern3)
  pat4 = re.compile(pattern4)
  pat5 = re.compile(pattern5)
  pat6 = re.compile(pattern6)
  pat7 = re.compile(pattern7)
  pat8 = re.compile(pattern8)
  pat9 = re.compile(pattern9)
  pattern = [pat1, pat2, pat3, pat4, pat5, pat6, pat7, pat8, pat9]
  f = open(filename, "r")
  line_index = 1
  space_pat = re.compile(r'}\s*')
  segment = []
  begin = []
  end = []
  flag = False
  rep_info = []
  head = False
  new_content = ''
  for line in f.readlines():
    print '\n\n'
    print 'origin:', line
    for (index, pat) in zip(range(0,9), pattern):
      res = pat.match(line)
      if res:
        res2 = list(res.groups())
        print 'origin res2:', res2
        for mat in new_list:
          if res2[1] == mat[0]:
            i = int(mat[1])
            flag_change = True
            print filename, line_index, res2[i + 2]
            for mat2 in new_list:
              if filename == mat2[2] and line_index > mat2[3] and line_index < mat2[4] and res2[i + 2] == mat2[5]:
                flag_change = False
            if flag_change:
              res2[i + 2] = '&' + res2[i + 2]
            flag = True
        print 'modified res2:', res2
    line2 = line
    if flag:
      line2 = '%s%s(' % (res2[0], res2[1])
      for ind2 in range(2, len(res2)):
        if ind2 == len(res2) - 1:
          line2 += res2[ind2] + ');\n'
        else:
          line2 += res2[ind2] + ','
      print 'line2:', line2
    new_content += line2
    flag = False
    line_index += 1
  f.close()
  f2 = open(filename, 'w')
  f2.write(new_content)
  f2.close()
print new_list
#print "test: "
#string = 'int func(int & a, int b) {'
#string = 'int keep_seismo(struct Receiver* , struct Wave, int it, int rank);'
#string = 'int Locate_Receiver(char * filename, struct Receiver *R, struct Coord *C, int rank);'
#for subpat in pattern:
#  pat = re.compile(subpat)
#  res = pat.match(string)
#  if res:
#    print res.groups()
