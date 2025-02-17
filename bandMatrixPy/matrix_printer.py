def print_1d(row,label=None, new_line=False, ret=False):
    s = []
    if label:
        s.append('%s = ' % label)
    s.append("| %s |\n" % "  ".join('%2.3f' % x for x in row))
    if new_line:
        s.append("\n")
    content = "".join(s)
    if ret:
        return content
    else:
        print(content,end='')
