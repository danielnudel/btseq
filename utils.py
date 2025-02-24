

LOG_FILE="log.txt"
CURRENT_LOGS_LEVEL=0


def log(text, logs_level=2):
    if CURRENT_LOGS_LEVEL > logs_level:
        return
    with open(LOG_FILE, 'a') as f:
        f.writelines(f"({logs_level}) {text}" + '\n')