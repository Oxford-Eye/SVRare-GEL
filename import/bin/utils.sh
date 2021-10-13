#!/usr/bin/bash
# retry until succeed. returns pid for kill afterward
# e.g. $0 ssh -i ~/.ssh/id_ed25519 -N -L 33061:127.0.0.1:3306 jyu@corp.gel.ac@phpgridzlogn001.int.corp.gel.ac & && mysqladmin --host='127.0.0.1' --port 33061 -u test -ptest processlist

function fail {
    echo $1 >&2
    exit 1
}

function retry {
    # retr \$1 until success, then run \$2
    local n=1
    local max=50
    local delay=5
    while true; do
        "$@" && break || {
            if [[ $n -lt $max ]]; then
                ((n++))
                echo "try failed. Attempt $n/$max:" >&2
                sleep $delay;
            else
                fail "The command has failed after $n attempts."
            fi
        }
    done
}

function is_mysql {
    port=$1; shift
    hostname=$1; shift

    psid=$(lsof -i :${port} | grep -m 1 'ssh' | tr -s ' ' | cut -d' ' -f2)
    if [ -z "$psid" ]; then
        echo 0
    else
        ps aux | grep ${psid} | grep -c ${hostname}
    fi
}


function find_port {
    starting_port=$1
    netstat --numeric --numeric-ports --listen | \
    awk -F'[ \t:]*' -v port_to_use=$starting_port '/[:]/ {
        a[i++] = $5
    }
    END {
        n = asort (a)
        for (i = 1; i <= n; i++)
            if (a[i] < port_to_use)
                continue
            else if (a[i] == port_to_use)
                port_to_use++
            else
                break;
        print port_to_use
    }'
}

function ssh_tunnel {
    local n=1
    local max=10
    local delay=5
    starting_port=$1; shift
    ssh_key=$1; shift
    host=$1; shift
    host_port=$1; shift
    ssh_socket=$1; shift
    
    port=$(find_port ${starting_port})
    while true; do
        ssh -i ${ssh_key} -M -S ${ssh_socket} -o ExitOnForwardFailure=yes -fnNT -L ${port}:127.0.0.1:${host_port} $host && echo ${port} && break || {
            if [[ $n -lt $max ]]; then
                ((n++))
                ((port++))
                port=$(find_port ${port})
                echo "try failed. Attempt $n/$max:" >&2
                sleep $delay;
            else
                fail "The command has failed after $n attempts."
            fi
        }
    done
}

function ssh_tunnel_old {
    local n=1
    local max=10
    local delay=5
    starting_port=$1; shift
    ssh_key=$1; shift
    host=$1; shift
    host_port=$1; shift
    ssh_socket=$1; shift
    
    port=$(find_port ${starting_port})
    ssh -i ${ssh_key} -M -S ${ssh_socket} -o ExitOnForwardFailure=yes -fnNT -L ${port}:127.0.0.1:${host_port} $host && echo ${port} || {
        mysql_running=$(is_mysql ${port} ${host})
        if [[ "$mysql_running" == 1 ]]; then
            echo ${port}
        else
            fail "The command has failed"
        fi
    }

}
