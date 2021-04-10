import uuid
import sys

UUID_BASE=uuid.UUID('00000000-0000-0000-0000-000000000000')

def convert(hashable: str) -> uuid.UUID:
    return uuid.uuid5(UUID_BASE, hashable)

def srr_string_to_sorted_string(srrs: str, 
                                splitchar: str=';',
                                joinchar: str=' ') -> str:
    srr_list = srrs.split(splitchar)
    srr_list = sorted(srr_list)
    return ' '.join(srr_list)

if __name__=='__main__':
    print(convert(srr_string_to_sorted_string(sys.argv[1])))
