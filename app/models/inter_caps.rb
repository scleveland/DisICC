class InterCap
  include DataMapper::Resource
  
  property :caps_id, Serial
  property :aasequence1_id, Integer, :required => true
  property :aasequence2_id, Integer, :required => true
  property :position_one, Integer, :required => true #this is the position of the aasequence1
  property :position_two, Integer, :required => true # this is the position of the partner in aasequence2
  property :mean_one, Float, :required => true
  property :mean_two, Float, :required => true
  property :correlation, Float, :required => true
  property :seq1_id, Integer, :required => true
  property :seq2_id, Integer, :required => true
  property :alignment1_id, Integer, :required => false
  property :alignment2_id, Integer, :required => false
  property :created_at, DateTime, :required => false
  
  def seq1
    Sequence.get(self.seq1_id)
  end
  
  def seq2
    Sequence.get(self.seq2_id)
  end
  
  def aa1
    AAsequence.get(aasequence1_id)
  end
  
  def aa2
    AAsequence.get(aasequence2_id)
  end
  
end