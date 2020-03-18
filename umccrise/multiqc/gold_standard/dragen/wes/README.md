```
https://aps2.platform.illumina.com/v1/workflows?Content-Type=application/json
{
     "name": "TumorNormalMultiple",
     "description": "A DRAGEN workflow for somatic variant calling on tumor and normal samples."
}
```

https://aps2.platform.illumina.com/v1/workflows/wfl.abb4ad43f7c24baf917c31720dfa5e20/versions
```
{
    "version": "muiltiqc_background",
    "definition": {
        "startat": "fastq-mapper",
        "states": {
            "fastq-mapper": {
                "type": "LaunchTask",
                "taskRun": {
                    "name": "{{Context.WorkflowRunId}}",
                    "description": "Step: fastqMapper",
                    "execution": {
                        "image": {
                            "name": "{{args.fastqMapper_imageName}}",
                            "tag": "{{args.fastqMapper_imageTag}}"
                        },
                        "command": "python3",
                        "args": [
                            "/mnt/map_fastqs.py",
                            "--fastq-csv",
                            "/data/FASTQ/fastq_list.csv",
                            "--tracking-sheet",
                            "/mnt/tracking-sheet.xlsx",
                            "--output-dir",
                            "/data/FASTQ/custom/"
                        ],
                        "inputs": [
                            {
                                "path": "/mnt/map_fastqs.py",
                                "url": "{{args.fastq_mapper_script}}",
                                "mode": "Download",
                                "type": "File"
                            },
                            {
                                "path": "/data/FASTQ/fastq_list.csv",
                                "url": "gds://{{args.fastq_list_volume}}/{{args.fastq_list_csv}}",
                                "mode": "Download",
                                "type": "File"
                            },
                            {
                                "path": "/mnt/tracking-sheet.xlsx",
                                "url": "gds://{{args.tracking_sheet_xlsx}}",
                                "mode": "Download",
                                "type": "File"
                            }
                        ],
                        "outputs": [
                            {
                                "path": "/data/FASTQ/custom/",
                                "url": "gds://{{args.fastq_list_volume}}/{{args.fastq_lists_outdir}}"
                            }
                        ],
                        "systemFiles": {
                            "url": "gds://{{args.log_outdir}}/{{Context.WorkflowRunId}}"
                        },
                        "environment": {
                            "resources": {
                                "type": "standardHiCpu",
                                "size": "medium"
                            }
                        },
                        "retryLimit": 0
                    }
                },
                "waitForCompletion": true,
                "Next": "retrieveFastqListCount"
            },
            "retrieveFastqListCount": {
                "Type": "CallApi",
                "Connection": "IAP",
                "Method": "GET",
                "Url": "/v1/folders",
                "UrlQuery": {
                    "volume.name": "{{args.fastq_list_volume}}",
                    "path": "/{{args.fastq_lists_outdir}}/SBJ*",
                    "include": "totalItemCount",
                    "pageSize": "1000"
                },
                "ResultFilterJs": "function filterResponse(responseContents) { var subjectIds = responseContents.items.map(function(item){ return item.name.split(\".\")[0];}); var result = { \"fastqListMaxIndex\": responseContents.itemCount - 1, \"ids\": subjectIds }; return result;}",
                "ResultPath": "$.subjects",
                "Next": "prepareDragenLoop"
            },
            "prepareDragenLoop": {
                "Type": "DefineResult",
                "Result": "{{$.subjects.fastqListMaxIndex}}",
                "ResultPath": "$.iterator",
                "Next": "iterateDragen"
            },
            "iterateDragen": {
                "Type": "Choice",
                "Choices": [
                    {
                        "Variable": "$.iterator",
                        "NumericLessThan": 0,
                        "Next": "preparePollLoop"
                    }
                ],
                "Default": "dragen"
            },
            "dragen": {
                "type": "LaunchTask",
                "resultPath": "$.taskRun",
                "taskRun": {
                    "name": "{{Context.WorkflowRunId}}",
                    "description": "Step: dragen (align, varcall, svcall)",
                    "execution": {
                        "image": {
                            "name": "{{args.dragen_imageName}}",
                            "tag": "{{args.dragen_imageTag}}"
                        },
                        "command": "bash",
                        "args": [
                            "-c",
                            "/opt/edico/bin/dragen --partial-reconfig DNA-MAPPER --ignore-version-check true; mkdir -p /ephemeral/ref/{{args.hg_assembly}}_{{args.dragen_imageTag}}; tar -C /ephemeral/ref/{{args.hg_assembly}}_{{args.dragen_imageTag}} -xvf /mount/index/{{args.hg_assembly}}/{{args.dragen_imageTag}}_ht.tar; /opt/edico/bin/dragen --lic-instance-id-location /opt/instance-identity -f --ref-dir /ephemeral/ref/{{args.hg_assembly}}_{{args.dragen_imageTag}} --fastq-list /data/LIST/fastqs_normal.csv --tumor-fastq-list /data/LIST/fastqs_tumor.csv --output-directory /output/results --output-file-prefix {{$.subjects.ids[{{$.iterator}}]}} --enable-map-align true --enable-map-align-output true --enable-variant-caller true --enable-duplicate-marking true --vc-enable-clustered-events-filter true --vc-enable-triallelic-filter true --enable-sv true"
                        ],
                        "inputs": [
                            {
                                "Path": "/mount/index/{{args.hg_assembly}}/{{args.dragen_imageTag}}_ht.tar",
                                "Url": "gds://umccr-refdata-dev/dragen/hsapiens/{{args.hg_assembly}}/{{args.dragen_imageTag}}_ht.tar"
                            },
                            {
                                "Path": "/data/LIST/fastqs_normal.csv",
                                "Url": "gds://{{args.fastq_list_volume}}/{{args.fastq_lists_outdir}}/{{$.subjects.ids[{{$.iterator}}]}}/{{$.subjects.ids[{{$.iterator}}]}}_normal.csv"
                            },
                            {
                                "Path": "/data/LIST/fastqs_tumor.csv",
                                "Url": "gds://{{args.fastq_list_volume}}/{{args.fastq_lists_outdir}}/{{$.subjects.ids[{{$.iterator}}]}}/{{$.subjects.ids[{{$.iterator}}]}}_tumor.csv"
                            },
                            {
                                "Path": "/mount/fastqs",
                                "Url": "gds://{{args.fastq_folder}}",
                                "type": "folder"
                            }
                        ],
                        "outputs": [
                            {
                                "path": "/output/results",
                                "url": "gds://{{args.result_dir}}/{{$.subjects.ids[{{$.iterator}}]}}"
                            }
                        ],
                        "systemFiles": {
                            "url": "gds://{{args.log_outdir}}/{{Context.WorkflowRunId}}"
                        },
                        "Environment": {
                            "resources": {
                                "Type": "fpga",
                                "Size": "medium"
                            }
                        },
                        "retryLimit": 0
                    }
                },
                "Next": "decrementIterator"
            },
            "decrementIterator": {
                "Type": "RunJavaScript",
                "Inputs": "{{$}}",
                "Script": "function action(input) {input.iterator = input.iterator - 1;if (typeof input.taskRuns === 'undefined') {input.taskRuns = [];}input.taskRuns.push(input.taskRun.TaskRunId);return input;}",
                "Next": "iterateDragen"
            },
            "preparePollLoop": {
                "Type": "DefineResult",
                "Result": "{{$.subjects.fastqListMaxIndex}}",
                "ResultPath": "$.iterator",
                "Next": "iteratePoll"
            },
            "iteratePoll": {
                "Type": "Choice",
                "Choices": [
                    {
                        "Variable": "$.iterator",
                        "NumericLessThan": 0,
                        "Next": "multiqc"
                    }
                ],
                "Default": "poll"
            },
            "poll": {
                "Type": "CallApiAndPoll",
                "ResultPath": "$.polling",
                "Method": "GET",
                "Connection": "IAP",
                "Url": "/v1/tasks/runs/{{$.taskRuns[{{$.iterator}}]}}",
                "PollingJs": "function processResponse(input) { var originalResults = JSON.parse(input.OriginalResults); var results = {}; if (input.ResponseMessage.IsSuccessStatusCode === false) { results.Action = 'Fail'; results.ErrorCause = input.MappedError; } else if (originalResults.status === 'Pending' || originalResults.status === 'Running') { results.Action = 'Poll'; results.NextPollingIntervalInSeconds = 1; } else { results.Action = 'Succeed'; } results.FinalResults = originalResults; return results; }",
                "Next": "decrementPollIterator"
            },
            "decrementPollIterator": {
                "Type": "RunJavaScript",
                "Inputs": "{{$}}",
                "Script": "function action(input){input.iterator=input.iterator-1; return input; }",
                "Next": "iteratePoll"
            },
            "multiqc": {
                "type": "LaunchTask",
                "taskRun": {
                    "name": "{{Context.WorkflowRunId}}",
                    "description": "Step: MultiQC InterOp + Dragen",
                    "execution": {
                        "image": {
                            "name": "{{args.multiqc_imageName}}",
                            "tag": "{{args.multiqc_imageTag}}"
                        },
                        "command": "multiqc",
                        "args": [
                            "/mnt/input",
                            "-o",
                            "/mnt/output",
                            "-m",
                            "dragen"
                        ],
                        "inputs": [
                            {
                                "path": "/mnt/input",
                                "url": "gds://{{args.result_dir}}",
                                "type": "folder"
                            }
                        ],
                        "outputs": [
                            {
                                "path": "/mnt/output",
                                "url": "gds://{{args.result_dir}}/multiqc"
                            }
                        ],
                        "systemFiles": {
                            "url": "gds://{{args.log_outdir}}/{{Context.WorkflowRunId}}"
                        },
                        "environment": {
                            "resources": {
                                "type": "standardHiCpu",
                                "size": "medium"
                            }
                        },
                        "retryLimit": 0
                    }
                },
                "waitForCompletion": true,
                "End": true
            }
        },
        "arguments": [
            {
                "name": "hg_assembly"
            },
            {
                "name": "fastqMapper_imageName",
                "value": "umccr/alpine_pandas",
                "Overridable": false
            },
            {
                "name": "fastqMapper_imageTag",
                "value": "1.0.1",
                "Overridable": false
            },
            {
                "name": "dragen_imageName",
                "value": "699120554104.dkr.ecr.us-east-1.amazonaws.com/public/dragen",
                "Overridable": false
            },
            {
                "name": "dragen_imageTag",
                "value": "3.5.7",
                "Overridable": false
            },
            {
                "name": "multiqc_imageName",
                "value": "umccr/multiqc_dragen",
                "Overridable": false
            },
            {
                "name": "multiqc_imageTag",
                "value": "1.2.1",
                "Overridable": false
            },
            {
                "name": "fastq_mapper_script",
                "value": "https://raw.githubusercontent.com/umccr/infrastructure/master/scripts/showcase/map_fastqs.py",
                "Overridable": false
            },
            {
                "name": "fastq_list_volume"
            },
            {
                "name": "fastq_lists_outdir"
            },
            {
                "name": "fastq_list_csv"
            },
            {
                "name": "tracking_sheet_xlsx"
            },
            {
                "name": "fastq_folder"
            },
            {
                "name": "result_dir"
            },
            {
                "name": "log_outdir"
            }
        ],
        "connections": [
            {
                "name": "IAP",
                "type": "IlluminaJwt",
                "host": "https://aps2.platform.illumina.com"
            }
        ]
    }
}
```

https://aps2.platform.illumina.com/v1/workflows/wfl.abb4ad43f7c24baf917c31720dfa5e20/versions/muiltiqc_background:launch
```
{
    "name": "muiltiqc_background_run7",
    "input": {
        "arguments": [
            {
                "name": "hg_assembly",
                "value": "hg38"
            },
            {
                "name": "fastq_list_volume",
                "value": "umccr-primary-data-dev"
            },
            {
                "name": "fastq_lists_outdir",
                "value": "vladsaveliev/multiqc_background/fastq_input_lists"
            },
            {
                "name": "fastq_list_csv",
                "value": "vladsaveliev/multiqc_background/inputs.csv"
            },
            {
                "name": "tracking_sheet_xlsx",
                "value": "umccr-primary-data-dev/vladsaveliev/multiqc_background/tracking_sheet.xlsx"
            },
            {
                "name": "fastq_folder",
                "value": "umccr-fastq-data-dev/PD/multiqc_controls"
            },
            {
                "name": "result_dir",
                "value": "umccr-primary-data-dev/vladsaveliev/multiqc_background/results"
            },
            {
                "name": "log_outdir",
                "value": "umccr-primary-data-dev/vladsaveliev/multiqc_background/logs"
            }
        ],
        "connections": [
            {
                "name": "IAP",
                "credentials": "{{ jwt }}"
            }
        ]
    }
}
```



