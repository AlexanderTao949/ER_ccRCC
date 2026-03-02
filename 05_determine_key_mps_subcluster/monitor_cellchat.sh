#!/bin/bash
# monitor_cellchat.sh

PID="1072227"  # 你的进程ID
LOGFILE="cellchat_monitor_$(date +%Y%m%d_%H%M%S).log"

echo "开始监控进程 $PID，日志保存至: $LOGFILE"
echo "时间戳,PID,CPU%,内存%,RSS(KB),VSZ(KB),状态" > "$LOGFILE"

while kill -0 $PID 2>/dev/null; do
    ps -p $PID -o pid,pcpu,pmem,rss,vsz,stat --no-headers | \
    awk -v ts="$(date '+%Y-%m-%d %H:%M:%S')" '{print ts","$1","$2","$3","$4","$5","$6}' >> "$LOGFILE"
    sleep 1
done

echo "进程 $PID 已结束，监控停止"
