#!/usr/bin/bash
#source
source /data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/script/setup
source /data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/TF-TWAS/script/util

# Global variables, including ones in the setup.sh
SCRIPT_PATH=$MODEL_PATH"/script"


# Main workflow
if [[ "$BASH_SOURCE" == "$0" ]] # to use the utilities directly from this script
then
    # arguments
    USAGE="try bash $0 help to learn more about how to use this program"
    use_function=${1:?$USAGE} # if not provide 1 arg, then print USAGE

    if [[ $use_function == "help" ]]
    then
        # call function
        display_helpmsg # print help page
    elif [[ $use_function == "setup" ]]
    then
       # use argument
       USAGE="bash $0 setup setup"
       setup_file=${2:?$USAGE}

       # call function
       setting_env $setup_file
       echo -e "DONE: setup"
    elif [[ $use_function == "run_example" ]] # this will call the util: run_elasticnetCV_TF_window_zero
    then
        # use argument
        USAGE="bash $0 run_example config prefix"
        config=${2:?$USAGE}
        prefix=${3:?$USAGE}

        # call function
        echo -e "::: call function: run_elasticnetCV_TF_window_zero :::"
        run_elasticnetCV_TF_window_zero $config $prefix
        rm ./runTWCV_TF_*
        echo -e "DONE: run_elasticnetCV_TF_window_zero"
    elif [[ $use_function == "run_elasticnetCV" ]]
    then
        # use arguments
        USAGE="bash $0 run_elasticnetCV config.R prefix.str"
        config=${2:?$USAGE}
        prefix=${3:?$USAGE}

        # call function
        echo -e "::: call function: run_elasticnetCV :::"
        run_elasticnetCV $config $prefix 
        rm ./runTWCV_TF_*
        echo -e "DONE: run_elasticnetCV"
    elif [[ $use_function == "run_elasticnetCV_TF_window_zero" ]]
    then
        # use arguments
        USAGE="bash $0 run_elasticnetCV_TF_window_zero config.R prefix.str"
        config=${2:?$USAGE}
        prefix=${3:?$USAGE}

        # call function
        echo -e "::: call function: run_elasticnetCV_TF_window_zero :::"
        run_elasticnetCV_TF_window_zero $config $prefix
        rm ./runTWCV_TF_*
        echo -e "DONE: run_elasticnet_TF_window_zero"
    elif [[ $use_function == "run_tf_background" ]]
    then
       # use arguments
       USAGE="bash $0 run_tf_background config.R prefix.str num_times"
       config=${2:?$USAGE}
       prefix=${3:?$USAGE}
       num=${4:?$USAGE}

       # call function
       echo -e "::: call function: run tf_background :::"

       out_dir=`grep "out_dir" $config | cut -d"=" -f2 | sed s/"'"//g`
       #out_dir_file="working_TW_TF_naive_model_exp_10-foldCV_elasticNet_alpha0.5_Illumina_"$prefix"_"$prefix".txt"
       out_dir_file=(working*$prefix.txt)
       output=$prefix".result.bg.txt"
       # iteration: num_times
       START=1
       END=$num
       pre_str=$prefix
       for i in $(eval echo "{$START..$END}")
       do
           echo -e "run background model the $i time"
           run_tf_background $config $pre_str"_"$i && cat $out_dir$out_dir_file >> $output
           echo -e "find result at $output"
           echo -e "DONE: run tf_background"
       done
    elif [[ $use_function == "run_tf_background_window_zero" ]]
    then
       # use arguments
       USAGE="bash $0 run_tf_background_window_zero config.R prefix.str num_times"
       config=${2:?$USAGE}
       prefix=${3:?$USAGE}
       num=${4:?$USAGE}

       # call function
       echo -e "::: call function: run tf_background_window_zero :::"

       out_dir=`grep "out_dir" $config | cut -d"=" -f2 | sed s/"'"//g`
       #out_dir_file="working_TW_TF_naive_model_exp_10-foldCV_elasticNet_alpha0.5_Illumina_"$prefix"_"$prefix".txt"
       out_dir_file=(working*$prefix.txt)
       output=$prefix".result.bg.window.zero.txt"
       # iteration: num_times
       START=1
       END=$num
       pre_str=$prefix
       for i in $(eval echo "{$START..$END}")
       do
           echo -e "run background model with window zero the $i time"
           run_tf_background_window_zero $config $pre_str"_"$i && cat $out_dir$out_dir_file >> $output
           echo -e "find result at $output"
           echo -e "DONE: run tf_background_window_size"
       done
    elif [[ $use_function == "generate_x_matrix" ]]
    then
        # use arguments
        USAGE="bash $0 generate_x_matrix expFile vcfFile outPath"
        expFile=${2:?$USAGE}
        vcfFile=${3:?$USAGE}
        outPath=${4:?$USAGE}
        check_file_extention $vcfFile "vcf"

        # call function
        echo -e "::: calling function: generate_x_matrix :::"
        outFiles=`generate_x_matrix $expFile $vcfFile $outPath`
        ordered_y=`echo $outFiles|cut -d, -f 1`
        ordered_x=`echo $outFiles|cut -d, -f 2`
        echo -e "find output files at $ordered_x"

        # checking sample order
        exp_sample_id=`check_file_cols $ordered_x 2 1`
        geno_sample_id=`check_file_cols $ordered_y 2 1`
        if [[ $exp_sample_id != $geno_sample_id ]]
        then
            echo -e "[Error] sample ordering are not consistent"
            echo "$exp_sample_id, $geno_sample_id"
            exit
        else
            rm $ordered_y
            echo -e "DONE: generate_x_matrix"
        fi
    elif [[ $use_function == "generate_y_matrix" ]]
    then
        # use arguments
        USAGE="bash $0 generate_y_matrix expFile vcfFile outPath"
        expFile=${2:?$USAGE}
        vcfFile=${3:?$USAGE}
        outPath=${4:?$USAGE}
        check_file_extention $vcfFile "vcf"

        # call function
        echo -e "::: calling function: generate_y_matrix :::"
        outFiles=`generate_y_matrix $expFile $vcfFile $outPath`
        ordered_y=`echo $outFiles|cut -d, -f 1`
        ordered_x=`echo $outFiles|cut -d, -f 2`
        echo -e "find output files at $ordered_y"

        # checking sample order
        exp_sample_id=`check_file_cols $ordered_x 2 1`
        geno_sample_id=`check_file_cols $ordered_y 2 1`
        if [[ $exp_sample_id != $geno_sample_id ]]
        then
            echo -e "[Error] sample ordering are not consistent"
            echo "$exp_sample_id, $geno_sample_id"
            rm $ordered_x
            exit
        else
            rm $ordered_x
            echo -e "DONE: generate_y_matrix"
        fi
    elif [[ $use_function == "generate_gene_annot" ]]
    then
        # use arguments
        USAGE="bash $0 generate_gene_annot posFile genelist prefix outPath"
        posFile=${2:?$USAGE}
        genelist=${3:?$USAGE}
        prefix=${4:?$USAGE}
        outPath=${5:?$USAGE}

        # call function
        echo -e "::: call function: generate_gene_annot :::"
        outFile=`generate_gene_annot $posFile $genelist $prefix $outPath`
        echo -e "find output file at $outFile"
        echo -e "DONE: generate gene_annot"
         




    # for case which the given utility name is not avilable
    else
    echo -e "InputError: need 1 argument, 0 is given"
    echo -e "See below for what's available in TF-TWAS:\n"
    display_helpmsg # print help page
    fi
fi
