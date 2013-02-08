task :read_email do
   gmail = Gmail.new("disiccapp","disorder for the win")
   gmail.inbox.emails.each do |email|
     subject = email.subject 
      if (subject.count "Disopred") > 0
        filepath= "temp_data/#{subject}"
         f = File.new(filepath, "w+")
         f.write(email.body)
         f.close
         #store_
         email.mark(:read)
         email.archive!
      else
        puts subject
        email.mark(:read) 
        email.archive!
      end
   end
end