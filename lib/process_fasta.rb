# Process an uploaded file for use with Resque
# 
# To use
#   Resue.enqueue(ProcessFastaSaveSequences, new_file)
class ProcessFasta
  @queue = :disicc
  
  def self.perform(new_file, fasta_hash, alignment_name, user_id)
    file = File.new(new_file, "r")
    counter = 0
    order_count = 0
    abrev_name = ""
    alignment_sequence = ""
    while (line = file.gets)
      if line.count(">") > 0 && counter > 0
        @sequence = Sequence.get(fasta_hash[abrev_name.strip])
        puts @sequence.abrev_name
        @alignment = Alignment.new(:seq_id => @sequence.id,
                         :alignment_name => alignment_name,
                         :align_order => order_count,
                         :alignment_sequence => alignment_sequence,
                         :fasta_title => abrev_name)
        if !@alignment.valid?
          puts @alignment.errors.inspect()
        end
        logger.debug "VALID"
        logger.debug { @alignment.errors.inspect }
        @alignment.save
        alignment_to_positions(@alignment)              
        #this is the sequene label
        abrev_name = line.gsub(">", "")

        order_count += 1
        alignment_sequence = ""
      elsif line.count(">") > 0
        abrev_name = line.gsub(">","")
        logger.debug { abrev_name }
        alignment_sequence =""
      elsif counter > 0
        alignment_sequence = alignment_sequence + line.lstrip.rstrip
      end
      counter = counter + 1
    end
     # puts fasta_hash[abrev_name]
     #  @sequence = Sequence.get(fasta_hash[abrev_name.strip])
     #  logger.debug "After"
     #  puts "After"  
     #  logger.debug "OHNO NONONONONO" 
     #  @alignment = Alignment.new(:seq_id => @sequence.seq_id,
     #                   :alignment_name => params[:alignment_name],
     #                   :align_order => order_count,
     #                   :alignment_sequence => alignment_sequence,
     #                   :fasta_title => abrev_name)
     #   logger.debug "SHAZZZZZAAAAAAAM"                 
     #  logger.debug { @alignment.valid?}
     #  if !@alignment.valid?
     #    puts @alignment.errors.inspect()
     #  end
     #  logger.debug "VALID"
     #  logger.debug { @alignment.errors.inspect }
     #  @alignment.save
     #  @alignment.alignment_to_positions              
     #  #this is the sequene label
     #  abrev_name = line.gsub(">", "")
     # 
     #  order_count += 1
     #  alignment_sequence = ""
    file.close

    # Message the user when action is complete.
    user = User.get(user_id)
    puts results
    puts "****** USER EMAIL #{user.email}  **************"
    puts DisMailer.email_user(user.email, "From DisICC:: Your New Alignment Job is complete for : #{@alignment.alignment_name}!")
  end
end