with open('data/QofT.dat', 'r') as file:
    lines = file.readlines()

processed_lines = []
empty_lines = []

for i, line in enumerate(lines):
    cleaned_line = line.replace('&', '').replace('/', '').strip()
    cleaned_line = cleaned_line.replace(', ,', ',')
    processed_lines.append(cleaned_line)
    if not cleaned_line:
        empty_lines.append(i)

header = processed_lines[0]

data_lines = []
for i, empty_line in enumerate(empty_lines):
    if i + 1 < len(empty_lines):
        mm = empty_lines[i] + 1
        nn = empty_lines[i+1]
        data_line = ''.join(processed_lines[mm:nn])
        data_lines.append(data_line)
    else:
        break

with open('data/QofT_formatted.dat', 'w') as file:
    file.write(header + '\n') 
    for data_line in data_lines:    
        file.write(data_line + '\n')