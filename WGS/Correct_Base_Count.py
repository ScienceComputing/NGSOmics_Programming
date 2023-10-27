def correct_count_base(input_base):
    base_count = {}

    k = 0
    while k < len(input_base):
        base = input_base[k]
        k += 1
        count = ''
        while k < len(input_base) and input_base[k].isdigit():
            count += input_base[k]
            k += 1
        count = int(count) if count else 1 
        base_count[base] = base_count.get(base, 0) + count

    compressed_string = ''.join([f'{base}:{count} ' for base, count in sorted(base_count.items())])
    return compressed_string

base = 'A57T20C3G10N2C50'
count_result = correct_count_base(base)
print(count_result) # A:57 C:53 G:10 N:2 T:20 
