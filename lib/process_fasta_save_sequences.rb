# Process an uploaded file for use with Resque
# 
# To use
#   Resue.enqueue(ProcessFastaSaveSequences, new_file)
class ProcessFastaSaveSequences
  @queue = :disicc
  
  def self.perform(new_file, user_id)
    File.open(new_file, "wb"){ |f| f.write(params['datafile'].read)}
    #sequence_list = Sequence.create_sequences_from_fasta(new_file, current_user.id)
    #sequences = Sequence.all(:seq_id => sequence_list.map{|s| s.seq_id})
    #sequences = []
    file = File.new(new_file, 'r')
    ff = Bio::FlatFile.new(Bio::FastaFormat, file)
    order_count = 0
    alignment = "s"
    ff.each_entry do |f|
      puts f.definition
      seq = Sequence.create_sequence_from_bioruby_fasta_entry(f, params[:seq_type],current_user.id)
      puts seq.abrev_name
      alignment = Alignment.create(:seq_id => seq.id,
                        :alignment_name => params[:alignment_name],
                        :align_order => order_count,
                        :alignment_sequence =>f.naseq)
      #alignment_to_positions(alignment) 
      alignment.alignment_to_positions  
      order_count += 1
    end
    # Run disorder
    alignment.run_and_store_disorder(thread_num=20)
    
    #Run CICP methods
    alignment.run_align_assess
    alignment.run_xdet_threaded(dir="temp_data", thread_num=4)
    alignment.import_xdet_new(thread_num=4)
    alignment.run_rate4site
    alignment.import_rate4site(thread_num=4)
    alignment.run_perl_caps
    alignment.import_caps_threaded(dir = "temp_data/#{self.alignment_name}/", thread_num=4)
    # Message the user when action is complete.
    user = User.get(user_id)
    puts results
    puts "****** USER EMAIL #{user.email}  **************"
    puts DisMailer.email_user(user.email, "From DisICC:: Your New Alignment and Analsis Job is complete for : #{alignment.alignment_name}!")
  end
end