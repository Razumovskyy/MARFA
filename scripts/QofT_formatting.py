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

print('empty_lines_indices: ', empty_lines)
print('len of processed lines: ', len(processed_lines))
input()

header = processed_lines[0]
print(header)
print('len of processed lines: ', len(processed_lines))
input()

data_lines = []
for i, empty_line in enumerate(empty_lines):
    if i + 1 < len(empty_lines):
        print(i, empty_lines[i], empty_lines[i+1])
        mm = empty_lines[i] + 1
        nn = empty_lines[i+1]
        print(mm, ':', nn)
        print('len of sliced processed line: ', len(processed_lines[mm:nn]))
        input()
        data_line = ''.join(processed_lines[mm:nn])
        data_lines.append(data_line)
    else:
        break

print(len(data_lines))

# Write to the output file
with open('data/QofT_formatted.dat', 'w') as file:
    file.write(header + '\n')  # Write the dimension line
    for data_line in data_lines:    
        file.write(data_line + '\n')  # Write the processed data line