#!/usr/bin/env nextflow

params.pwd = "$PWD"
params.input = ""$params.pwd"/data"

process createFolders {

    script:
    """
    mkdir -p $params.pwd/log
    mkdir -p $params.pwd/tmp
    mkdir -p $params.pwd/output
    """

}

process startServer {

    script:
    """
    python -m client_server.live_server
    """

}

process startWatcher {

   output:
   stdout into result2

   script:
   """
   python -m watcher.watcher "$params.watched_dir"
   """
}


