class AAsequence 
  include DataMapper::Resource
  
  property :id, Serial, :field => 'AAsequence_id'
  property :seq_id, Integer, :required => true
  property :amino_acid, String, :length=> 1,  :required => true
  property :original_position, Integer, :required => false
  property :disorder_consensus, Float, :required => false, :default => 0.0
  property :contact_consensus, Float, :required => false, :default => 0.0
  property :contact_positive_consensus, Integer, :required => false, :default => 0
  property :deleted_at, ParanoidDateTime
  
  alias :AAsequence_id :id
  belongs_to :sequence, 'Sequence', :child_key => [:seq_id], :parent_key => [:seq_id]
  #has n, :disorder, 'Disorder', :parent_key=>[:disorder_id]
  #has n, :disorder_values, 'DisorderValue', :child_key => [:aasequence_id], :parent_key =>[:id]
  #has 1, :xdet, 'Xdet', :child_key => :aasequence_id
  #has 1, :conseq, 'Conseq', :child_key => :aasequence_id
  
  
  def calculate_intra_consensus_value
    count=0
    if IntraResidueContact.all(:seq_id => self.seq_id, :first_residue=>self.original_position) || IntraResidueContact.all(:seq_id => self.seq_id, :second_residue=>self.original_position)
      count +=1
    end
    if !Conseq.first(:aasequence_id => self.AAsequence_id).nil? && Conseq.first(:aasequence_id => self.AAsequence_id).color > 4
        count +=1
    end
    if !Xdet.first(:aasequence_id => self.AAsequence_id).nil? && (Xdet.first(:aasequence_id => self.AAsequence_id).correlation > 0.0 || Xdet.first(:aasequence_id => self.AAsequence_id).correlation == -2)
        count +=1
    end
    if !NewCap.first(:seq_id=> self.seq_id, :position_one => self.original_position).nil? || !NewCap.first(:seq_id=> self.seq_id, :position_two => self.original_position).nil?
      count +=1
    end
    self.contact_consensus = count /4
    puts self.contact_consensus
    self.save
  end
  
  def calculate_intra_consensus_value_special
    count=0
    if !Conseq.first(:aasequence_id => self.AAsequence_id).nil? && Conseq.first(:aasequence_id => self.AAsequence_id).color > 4
        count +=1
    end
    if !Xdet.first(:aasequence_id => self.AAsequence_id).nil? && (Xdet.first(:aasequence_id => self.AAsequence_id).correlation > 0.0 || Xdet.first(:aasequence_id => self.AAsequence_id).correlation == -2)
        count +=1
    end
    if !NewCap.first(:seq_id=> self.seq_id, :position_one => self.original_position).nil? || !NewCap.first(:seq_id=> self.seq_id, :position_two => self.original_position).nil?
      count +=1
    end
    self.contact_consensus = count/3
    self.save
    puts self.contact_consensus
  end
  
end
